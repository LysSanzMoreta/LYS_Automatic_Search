# LYS_Automatic_Search

Please feel free to submit an issue if you encounter any :) and a leave a like if these scripts were useful tou you. Keep in mind that these scripts are under low maintenance.

**Citation**: If you use this scripts make sure to cite : https://link.springer.com/article/10.1007/s11692-020-09507-9 <br/>


**Summary**
This script can perform the following actions: 

i) Mines the RCSB database [10] using the blast alignment tool to find the best matching homologous sequences 
ii) Fetch their domain positions by using Prosites [3,8,9] 
iii) Parse the output of PAML extracting the positional information of fast-evolving sites and transform them into the coordinate system of the protein structure iv) Output a file per gene with the positions correlations to its homologous sequence. This output is used as an input in the GUi to generate the graphical assesment. <br/>


WARNING: TKINTER STRUGGLES WITH HIGH DEFINITION SCREENS, not resolved yet by the tkinter people. So the GUI will look strange
WARNING: I have simplidied the GUI from the original paper, since it is not necessary anymore. It is still under construction for improvements.

Expected result: ![alt text](https://github.com/artistworking/LYS_Automatic_Search/blob/master/MKKS_2.png)

**Python environment**: 

LYS_Automatics_Search/environment.yaml

**Example**: Run LYS_autosearch_example.py together with the data under -TestSequences- for pipeline testing

Required input files:

- proteins: Path to fasta file with the gene sequences (i.e Test5.fasta)
- codeml-output: Path to folder with the codeml outputs (i.e TestSequences)

Optional input files:

- blast-results: Once the autosearch has been run at least once, you will obtain a file named Full_Blast_results_against_PDB.tsv, you can re-use it here
- pdb-files: The PDB files will be store here and can also be re-used
- pdb-results: Once the autosearch has been run at least once, you will obtain a file named Full_Blast_results_against_PDB_Filtered.tsv, you can re-use it here
 
Find the other arguments explained in the script


