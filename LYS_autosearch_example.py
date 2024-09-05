import os,sys,ntpath
import argparse
import datetime
from argparse import RawTextHelpFormatter
local_repository=True
script_dir = os.path.dirname(os.path.abspath(__file__))
if local_repository:
     sys.path.insert(1, "{}/lysautosearch/src".format(script_dir))
     import lysautosearch
else:#pip installed module
     import lysautosearch
import lysautosearch as LASutils
now = datetime.datetime.now()
def main(args):
    results_dir = "{}/RESULTS_LysAutoSearch_{}_{}".format(script_dir,args.results_name,now.strftime("%Y_%m_%d_%Hh%Mmin%Ss%fms"))
    args.__dict__["results_dir"] = results_dir
    basepath = ntpath.basename(results_dir)
    LASutils.folders(basepath, script_dir)
    if not os.path.exists(args.pdb_files):
        LASutils.folders("PDB_files", script_dir,overwrite=False) #only overwrite=True if we want to re-download everything
        args.__dict__["pdb_files"] = "{}/{}".format(script_dir,"PDB_files")
    lysautosearch.run(args,parser,results_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="LYS auto search args",formatter_class=RawTextHelpFormatter)
    # COMPULSORY files: Give error if not present
    parser.add_argument('--results-name',
                        help='Run/Results name',default="Test")
    parser.add_argument('--proteins',
                        help='FULL Path to file with all positively selected peptide sequences',default="TestSequences/Test5.fasta")
    parser.add_argument('--codeml-output',
                        help='FULL path to tab separated text file with the fasta ids and the paths to codeml output files (M8)',default="TestSequences/Gene_paths.txt")
    # OPTIONAL arguments:
    parser.add_argument('--blast-results',
                        help='Path to previously run Blast results dataframe found under the name of -Full_Blast_results_against_PDB.tsv-',
                        default="/home/lys/Dropbox/MasterRepositories/LYS_Automatic_Search/RESULTS_LysAutoSearch_Test_2024_09_05_19h24min52s188505ms/Full_Blast_results_against_PDB.tsv")

    parser.add_argument('--pdb-files',
                        help='Path to folder with PDB files, otherwise a new folder will created and new download run will take place',
                        default="/home/lys/Dropbox/MasterRepositories/LYS_Automatic_Search/PDB_files")
    parser.add_argument('--pdb-results',
                        help='Path to file containing the matrix with the PDB files from the homologous sequences under the name of -Full_Blast_results_against_PDB_Filtered.tsv-',
                        #default="",
                        default="/home/lys/Dropbox/MasterRepositories/LYS_Automatic_Search/RESULTS_LysAutoSearch_Test_2024_09_05_19h24min52s188505ms/Full_Blast_results_against_PDB_Filtered.tsv"
                        )


    parser.add_argument('--fasta-format',
                        help='Sequence or Multiple Alignment File Format (default fasta), do not write it with quotes',
                        default='fasta')
    parser.add_argument('--number-homologous',
                        help='Maximum number of homologous sequences to keep',
                        default=3)
    parser.add_argument("--probability",
                        help='Choice of level of posterior probability on the sites, 95% or 99% from M8 Out file',
                        default=int(
                             99))  # prob : Posterior probability percentage of the positive sites , 99 or 95, default= 99
    parser.add_argument("--missing_data",
                        help='Decide if the missing data("N") should be removed from the nucleotide sequence, default = yes, recommended for the alignment',
                        default='yes')  # missing_data : remove or not missing data from the sequence alignment: yes or no, default=yes
    parser.add_argument("--print_alignment",
                        help='Choose to visualize the PDB file sequence aligned with the gene, default = no',
                        default='no')
    parser.add_argument("--number_homologous",
                        help='Select the top n homologous proteins to be chosen to perform the positions dataframes ordered by resolution',
                        default=int(3))

    args = parser.parse_args()
    args.__dict__["script_dir"] = script_dir


    main(args)

