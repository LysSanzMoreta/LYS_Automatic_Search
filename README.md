# LYS_Automatic_Search




Citation: If you use this scripts make sure to cite : https://link.springer.com/article/10.1007/s11692-020-09507-9 <br/>


Scripts that i) Mine the RCSB database [10] using the blast alignment tool to find the best matching homologous sequences ii) Fetch their domain positions by using Prosites [3,8,9] iii) Parse the output of PAML extracting the positional information of fast-evolving sites and transform them into the coordinate system of the protein structure iv) Output a file per gene with the positions correlations to its homologous sequence. This output is used as an input in the GUi to generate the graphical assesment. The basic tutorial on these scripts is available at https://youtu.be/ZUxUHWfQ9kw<br/>

Requires: Python3, Pandas, Numpy, Tkinter, PyMOL (I am afraid currently a license is needed to use the latest pyMOL that works in python 3) ---> Use anaconda 3 to install easily

WARNING: TKINTER STRUGGLES WITH HIGH DEFINITION SCREENS, not resolved yet by the tkinter people. So the GUI will look strange


![#f03c15](https://placehold.it/15/f03c15/000000?text=+) FIRST STEP:<br/>
-LYS_PDB_Search.py:<br/>
    ```python3 LYS_PDB_Search.py --Genes Proteins/Trial_sequences.fasta --Codeml Gene_and_path_to_codeml_result.fasta```
    <br/>
    a) INPUT FILES:<br/>
        1) File with coding sequences to analyze: Whatever format, preferably, simple identifiers without spaces (to avoid errors) <br/> 
        2) Tab separated file : Column 1: gene identifier (same as in file a)) and Column 2: Path to the codeml output file containing the + sites <br/>
    b)OUTPUT FOLDERS:<br/>
        1)PDB_files: Contains the automatically downloaded PDB files (according only to the filtered results, saving memory and time)<br/>
        2)Positions_Dataframes: Contains the dataframes with the corresponding positions from the gene in the PDB sequence---> Offers the most important information: Search for 'Selected_and_Domain' (the most interesting positions where the +site is in the domain): <br/> 
        
3)Lys_Pymol_Images: Empty for now <br/>

![#f03c15](https://placehold.it/15/f03c15/000000?text=+) SECOND STEP:<br/>
-LYS_PyMOL_input_Dataframe_GUI.py: <br/> If using in server and pymol python package version 2.2 installed via anaconda3 does not work download files outside the server and run the script with pymol python package version 1.8. Current versions of Pymol work too but you require a license <br/>
```python3 LYS_PyMOL_input_Dataframe_GUI.py```  <br/>
Use this command sequence as examples to extract the Position dataframes and their corresponding PDB files to separated folders (again, in case pymol does not work properly in the server):<br/>
        
        ```bash
        #In the LYS_PDB_Search.py directory
        mkdir Dataframes_Selected_Domains
        cd Positions_Dataframe
        cp `grep -lir 'Selected_and_Domain' .` ../Dataframes_Selected_Domains
        cd ..
        ls Dataframes_Selected_Domains >> Filenames.txt
        sed 's/[^_]*_\([^_]*\)_.*/\1/g' Filenames.txt > Filenames2.txt
        mkdir PDB_files_Selected_and_Domains
        cd /home/lys/HUMMING/PDB_files
        cat ../Filenames2.txt |xargs cp -t ../PDB_files_Selected_and_Domains/
        ```
  a)INPUT FILE:<br/>
        Introduce in the GUI the path to the dataframe that the user wants to plot from the folder Positions_Dataframe. Adjust conveniently <br/>
        ATTENTION: If you executed and downloaded from a server the folders Dataframes_Selected_Domains and PDB_files_Selected_and_Domains, and you are working locally. Remember to rename the second folder to PDB_files and always work in the directory where the Dataframes_Selected_Domains and PDB_files are!! As well re-create the Lys_Pymol_Images folder , mimicking the one in the server (It will not be able to find the PDB files otherwise or where to put the images) <br/>
  b)OUTPUT FILE: <br/>
        Lys_Pymol_Images: The images will be directed to this folder<br/>
![alt text](https://github.com/artistworking/LYS_Automatic_Search/blob/master/MKKS_2.png)
        
    Examples of commands:

    # Test multiple sequences
    python3 LYS_PDB_Search.py --Genes Proteins/Trial_sequences.fasta --Codeml Gene_and_path_to_codeml_result.fasta
    # Test 1 sequence with a PDB file of personal choice (not via automatic search), the domains detection CAN be automatic (Prosites) or not
    python3 LYS_PyMOL_Prosites --PDB TLR7.pdb --Gene GENEaln/New_IDS/ENSTGUT00000008587.map2finch.fasta --M8 TLR7_ENSTGUT00000008587.map2finch.fasta/M8/out --prosite yes --sequence_number 1 --missing_data no
    #GUIs to test individually each protein
    #a) with dataframes obtained from LYS_PDB_Search.py
    python3 LYS_PyMOL_input_Dataframe_GUI.py
    #b) GUI equivalent to LYS_PyMOL_Prosites.py: Test 1 sequence with a PDB file of personal choice (not via automatic search), the domains detection CAN be automatic (Prosites) or not (if NOT, you need to download and input manually a fasta file with the domain sequences (caution cause they are collapsed into one, check them one by one) and it also needs the full PDB sequence as referred in this old version of the same script video tutorial: https://www.youtube.com/watch?v=v3MawxHLGSo , but the idea is the same)
    python3 LYS_PyMOL_GUI_Prosites.py
    
Additional Comments: <br/>

If the number of input sequences and the number of output dataframes is not the same is due to Prosites not been able to find functional domains on the structure. It is relatively common, therefore, simply use the additional scripts (LYS_PyMOL_GUI_Prosites.py or LYS_PyMOL_Prosites.py) to download the functional domains sequences from somewhere else (eg. Uniprot) and use the recommended PDB file from Full_Blast_results_against_PDB_Filtered.tsv. https://www.youtube.com/watch?v=v3MawxHLGSo
        
        
 Have fun!!! :D
