import requests
import time
import argparse
import sys
import os
import json
import csv
from pathlib import Path

BASE = "https://wwwdev.ebi.ac.uk/pdbe/pdbe-kb/pisa/api"


def submit_pdb(pdb_path):
    url = f"{BASE}/submit"

    with open(pdb_path, "rb") as f:
        files = {"file": f}
        print(f"[DEBUG] POST {url} :: {pdb_path.name}")
        resp = requests.post(url, files=files, timeout=60)

    if not resp.ok:
        print(resp.status_code)
        print(resp.text)
        resp.raise_for_status()

    data = resp.json()
    if "job_id" not in data:
        raise RuntimeError(f"job_id not found in response: {data}")

    return data["job_id"]


def fetch_interface_summary(job_id):
    url = f"{BASE}/results/interface_summary/{job_id}"
    print(f"[DEBUG] GET {url}")
    resp = requests.get(url, timeout=60)
    return resp


def wait_for_results(job_id, interval=5, timeout=600):
    start = time.time()

    while True:
        resp = fetch_interface_summary(job_id)

        if resp.status_code == 200:
            return resp.json()

        if resp.status_code in (202, 204, 404):
            if time.time() - start > timeout:
                raise TimeoutError("PISA results timeout")
            time.sleep(interval)
            continue

        print(resp.status_code)
        print(resp.text)
        resp.raise_for_status()


def _normalize_record(record):
    area_keys = ["int_area", "interface_area", "area"]
    solv_keys = ["int_solv_energy", "solvation_energy", "solv_energy", "delta_g"]
    pvalue_keys = ["pvalue", "p_value", "p-val", "p_val"]

    area = next((record.get(k) for k in area_keys if k in record), None)
    solv = next((record.get(k) for k in solv_keys if k in record), None)
    pvalue = next((record.get(k) for k in pvalue_keys if k in record), None)

    return {
        "int_area": area,
        "int_solv_energy": solv,
        "pvalue": pvalue,
    }


def _looks_like_interface_record(obj):
    if not isinstance(obj, dict):
        return False

    keys = set(obj.keys())
    return any(
        k in keys
        for k in [
            "int_area",
            "int_solv_energy",
            "pvalue",
            "interface_area",
            "solvation_energy",
            "delta_g",
        ]
    )


def _walk_for_interface_records(obj, found):
    if isinstance(obj, dict):
        if _looks_like_interface_record(obj):
            found.append(_normalize_record(obj))
        for value in obj.values():
            _walk_for_interface_records(value, found)
    elif isinstance(obj, list):
        for item in obj:
            _walk_for_interface_records(item, found)


def extract_metrics(data):
    results = []

    if isinstance(data, dict) and isinstance(data.get("interfaces"), list):
        for interface in data["interfaces"]:
            if isinstance(interface, dict):
                results.append(_normalize_record(interface))

    if not results:
        _walk_for_interface_records(data, results)

    deduped = []
    seen = set()
    for item in results:
        key = (
            item.get("int_area"),
            item.get("int_solv_energy"),
            item.get("pvalue"),
        )
        if key not in seen:
            seen.add(key)
            deduped.append(item)

    return deduped


def save_debug_json(data, job_id, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, f"{job_id}_interface_summary_raw.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)
    print(f"[DEBUG] Saved raw response to {path}")


def find_pdb_files(input_dir):
    input_path = Path(input_dir)
    if not input_path.exists():
        raise FileNotFoundError(f"Directory not found: {input_dir}")
    if not input_path.is_dir():
        raise NotADirectoryError(f"Not a directory: {input_dir}")

    pdb_files = sorted(input_path.glob("*.pdb"))
    return pdb_files


def run_single_pdb(pdb_path, output_dir):
    job_id = submit_pdb(pdb_path)
    print(f"[INFO] Job ID for {pdb_path.name}: {job_id}")

    data = wait_for_results(job_id, interval=5, timeout=600)
    save_debug_json(data, job_id, output_dir)

    results = extract_metrics(data)
    return job_id, results


def write_csv(rows, csv_path):
    fieldnames = [
        "pdb_file",
        "job_id",
        "interface_index",
        "int_area",
        "int_solv_energy",
        "pvalue",
    ]

    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description="Run PDBe-KB PISA for all PDB files in a directory")
    parser.add_argument("--input", required=True, help="Input directory containing PDB files")
    parser.add_argument(
        "--output-dir",
        default="pisa_output",
        help="Directory for raw JSON outputs and final CSV",
    )
    args = parser.parse_args()

    try:
        pdb_files = find_pdb_files(args.input)
    except Exception as exc:
        print(f"[ERROR] {exc}")
        sys.exit(1)

    if not pdb_files:
        print(f"[ERROR] No .pdb files found in directory: {args.input}")
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)
    csv_path = os.path.join(args.output_dir, "PISA_interface.csv")
    all_rows = []

    print(f"[INFO] Found {len(pdb_files)} PDB file(s)")

    for idx, pdb_path in enumerate(pdb_files, start=1):
        print("=" * 80)
        print(f"[INFO] Processing {idx}/{len(pdb_files)}: {pdb_path.name}")

        try:
            job_id, results = run_single_pdb(pdb_path, args.output_dir)

            if not results:
                print(f"[WARN] No interface records parsed for {pdb_path.name}")
                all_rows.append(
                    {
                        "pdb_file": pdb_path.name,
                        "job_id": job_id,
                        "interface_index": "",
                        "int_area": "",
                        "int_solv_energy": "",
                        "pvalue": "",
                    }
                )
                continue

            for interface_index, result in enumerate(results, start=1):
                print(
                    f"[RESULT] {pdb_path.name} | Interface {interface_index}: "
                    f"int_area={result['int_area']}, "
                    f"int_solv_energy={result['int_solv_energy']}, "
                    f"pvalue={result['pvalue']}"
                )

                all_rows.append(
                    {
                        "pdb_file": pdb_path.name,
                        "job_id": job_id,
                        "interface_index": interface_index,
                        "int_area": result["int_area"],
                        "int_solv_energy": result["int_solv_energy"],
                        "pvalue": result["pvalue"],
                    }
                )

        except Exception as exc:
            print(f"[ERROR] Failed for {pdb_path.name}: {exc}")
            all_rows.append(
                {
                    "pdb_file": pdb_path.name,
                    "job_id": "",
                    "interface_index": "",
                    "int_area": "",
                    "int_solv_energy": "",
                    "pvalue": "",
                }
            )

    write_csv(all_rows, csv_path)
    print("=" * 80)
    print(f"[INFO] Final CSV saved to: {csv_path}")


if __name__ == "__main__":
    main()
