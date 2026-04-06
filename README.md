<p align="center">
<img src="./images/mint_choco_latte.png" alt="cats" width="75%"/>
</p>

## Table of Contents
- [af3filter](#af3filter)
- [mpnn2afserver](#mpnn2afserver)
- [pdb2af3server](#pdb2af3server)
- [Run PISA interface](#Run-PISA-interface)
- [pdb2png](#pdb2png)
- [af3db_rename](#af3db_rename)
  
</br>

**Installation:**

Linux
```bash
git clone https://github.com/scottcsh/protdesign
```
Windows
<pre>
Code > Download ZIP
</pre>

# af3filter
Filter AlphaFold 3 Server Outputs

</br>

**Requirements:**

jq in Centos7
```bash
yum install jq -y
```

</br>

**Usage:**


1. Run & Download alphafold results from https://alphafoldserver.com/ in zip file

	(or run alphafold 3 in the local workstation)

2. Unzip the results, run the script as:
```bash
./AF3_filter.sh --dir <json directory> --out <output filename> --iptm <threshold> --ptm <threshold> --score <threshold> --pae <threshold> --chain_id <A-E> --max_output <N>
```

</br>

**Options:**
<pre>
  --dir			Root directory containing JSON files (required)
  
  --out			Output CSV filename (optional, default: results.csv)
  
  --iptm		Minimum iptm value (optional)
  
  --ptm			Minimum ptm value (optional)
  
  --score		Minimum ranking_score value (optional)
  
  --pae			Maximum allowed chain_pair_pae_min value (optional)
  
  --chain_id    Chain ID (A–E) to select chain_pair_pae_min value (required for --pae option)
  
  --max_output	Limit output to top N entries sorted by average(iptm, ptm, score)
  
  --help		Show help message
</pre>

</br>

**Example:**
```bash
./AF3_filter.sh --dir folds_2026_01_20_02_08/ --iptm 0.7 --ptm 0.7 --score 0.8 --pae 1.8 --chain_id A --max_output 50
```

</br>

# mpnn2afserver
Process ProteinMPNN fasta files into AFserver_input.json

</br>

**Usage:**

```bash
./mpmm2afserver.sh --dir <fasta directory> --fa <target protein fasta> --max_job <number of jobs per json>
```
</br>

**Example:**
```bash
./mpnn2afserver.sh --dir seqs/ --fa 1sy6.fa --max_job 100
```

</br>

# pdb2af3server
Process pdb files into AFserver_input.json

</br>

**Usage:**

```bash
./pdb2af3server.sh --dir <pdb directory> --target <target protein pdb> --chain <chain id>
```
</br>

**Example:**

```bash
./pdb2af3server.sh --dir seqs/ --target 1sy6.pdb --chain A
```
</br>

# Run PISA interface
Run PISA on pdb files

</br>

**Usage:**

```bash
python run_pisa_interface.py --dir <pdb directory>
```
</br>

**Example:**

```bash
./pdb2af3server.sh --dir folder_with_pdbs/
```
</br>


# pdb2png

</br>

**Requirements:**

PyMOL

</br>

**Usage:**

1. Edit **pdb2png_PyMOL.py** with notepad
   
   <pre>
   pdb_path = r"C:\dir\input\*.pdb"
   output_dir = r"C:\dir\folder_for_output_thumbnails"
   </pre>
2. Run **PyMOL**
3. File > Run Script > **pdb2png_PyMOL.py**

</br>

# af3db_rename
Rename seqkit split2 results of sharded alphafold3 genetic databases
</br>
Reference: (https://github.com/google-deepmind/alphafold3/blob/main/docs/performance.md#sharded-genetic-databases)
</br>

**Usage:**

```bash
./AF3DB_rename.sh --dir <dir_for_af3db_split.fasta>
```
</br>

**Example:**


**Input**
<pre>
	bfd-first_non_consensus_sequences_shuffled.part_001.fasta
	bfd-first_non_consensus_sequences_shuffled.part_002.fasta
							...
	bfd-first_non_consensus_sequences_shuffled.part_064.fasta
</pre>

</br>

**Output**
<pre>
	bfd-first_non_consensus_sequences_shuffled.fasta-00000-of-00064
	bfd-first_non_consensus_sequences_shuffled.fasta-00001-of-00064
							...
	bfd-first_non_consensus_sequences_shuffled.fasta-00063-of-00064
</pre>


</br>
</br>

[Return to top](#Table-of-Contents)
