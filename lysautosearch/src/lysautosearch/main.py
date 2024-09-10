import os.path
import pandas as pd
import lysautosearch
import lysautosearch.utils as LASutils
import lysautosearch.blast_tools as LASblast
from Bio import SeqRecord, SeqIO

def args_checker(args,parser):
    # Check all the necessary files are present and the right arguments are given:
    file_formats = ['fasta', 'phylip_sequential', 'clustal', 'embl', 'genebank', 'gb', 'abi', 'ace', 'fastq-sanger',
                    'fastq', 'fastq-solexa', 'fastq-illumina', 'ig', 'imgt', 'pdb-seqres', 'pdb-atom', 'phd', 'phylip',
                    'pir', 'seqxml', 'sff', 'stockholm', 'swiss', 'tab', 'qual', 'uniprot-xml']
    if not args.proteins or not args.codeml_output: #TODO: Refactor Genes to proteins
        parser.error('Missing one of the 2 required arguments: Genes sequences / Dataframewith codeml paths and gene name')
        # Genes = "/home/lys/Downloads/LYS_Automatic_Search-master/TestSequences/Test5.fasta"
        # codeml_output = "/home/lys/Downloads/LYS_Automatic_Search-master/TestSequences/Fasta_paths.txt"
    elif args.probability not in [95, 99]:
        parser.error('Probability values can only be 95  or 99')
    elif args.missing_data not in ['no', 'yes']:
        parser.error('Choose to keep missing data: yes or to remove it: no')
    elif args.fasta_format not in file_formats:
        parser.error('Invalid Gene File Format.Check available SeqIO biopython file formats ')
    else:
        print('Building dataframe...')


def start_blast(args,results_dir):
    # Highlight: Refactored List_protein_sequences to protein_sequences
    protein_sequences, protein_id_names = LASutils.fasta_to_sequences(args.proteins,
                                                                      args.fasta_format)  # highlight: refactored Genes to args.proteins
    protein_sequences = list(map(lambda seq: (LASutils.translate_sequence(seq,remove_missing_data=False),LASutils.translate_sequence(seq,remove_missing_data=True)), protein_sequences))

    protein_sequences = list(zip(*protein_sequences))
    protein_sequences_missing_data = protein_sequences[0]
    protein_sequences_no_missing_data = protein_sequences[1]
    # Highlight: Generate the full and filtered summary dataframe with the top suggested homologous proteins from PDB
    blast_df = LASblast.blast_dataframe(protein_sequences_missing_data,protein_sequences_no_missing_data , protein_id_names, results_dir)


    return blast_df

#def Wrapper_of_all_functions(pdb_sequence,Gene,Chains,M8,List_Domains,Format,prob,missing_data,Residues_ID,basepath,print_alignment,Gene_name):
def build_PyMOL_dataframe(pdb_sequence, sequences_dict, gene_name, chains, M8_file, domains_idx, residues_ID, pdb_file_name, args):
    '''Calls the function that builds the dataframe of positions that will be read by PyMOL according to according to the optional arguments
    given by the user'''

    M8_file = os.path.join(args.codeml_output,M8_file)

    positives_sites =LASutils.extract_positive_sites(M8_file,args.probability)
    prot_missing_data= sequences_dict["prot_missing_data"]
    prot_no_missing_data= sequences_dict["prot_no_missing_data"]
    #Checking if the user wants to perform the alignment with or without missing data in the gene
    if args.missing_data == 'no': ###The sequence 'Gene' will be translated at this point no matter what
        #TODO: refactor Clean_positions to mapped_positions
        mapped_positions = LASutils.corresponding_positions_missing_notmissing_data(prot_missing_data,prot_no_missing_data) #Find the equivalent positions among the protein with and without missing data
        #pymol_dataframe =Corresponding_Coordinates_and_labels_PDB_Gene(pdb_sequence,prot_no_missing_data, positives_sites,domains_idx,residues_ID,basepath,args.print_alignment,gene_name, mapped_positions)
        pymol_dataframe =LASutils.corresponding_coordinates_and_labels_PDB_gene(pdb_sequence, pdb_file_name,prot_no_missing_data, positives_sites,domains_idx,residues_ID,args.proteins,args,gene_name, mapped_positions)

    else:
        #pymol_dataframe =Corresponding_Coordinates_and_labels_PDB_Gene(pdb_sequence,prot_missing_data,List_of_Positive_Positions,List_Domains,Residues_ID,basepath,print_alignment,Gene_name)

        pymol_dataframe = LASutils.corresponding_coordinates_and_labels_PDB_gene(pdb_sequence, pdb_file_name, prot_missing_data,
                                                                                 positives_sites, domains_idx,
                                                                                 residues_ID, args.proteins, args,
                                                                                 gene_name)
        #Corresponding_Coordinates_and_labels_PDB_Gene(PDB,Gene,List_of_Positive_Positions,functional_list,Residues_ID,basepath,print_alignment,Clean_positions = None):
    return pymol_dataframe

def run(args,parser,results_dir): #transferring pdb_search here
    """Blast each gene against the RSCB database and retrieve the PDB file from homologous proteins.
     Call Pymol and color the homologous proteins by the positively selected residues indicated by PAML
     """
    #Read the input fasta file and divide the gene ids and the sequences
    args_checker(args,parser)

    if not args.blast_results:
        blast_df = start_blast(args,results_dir)
    else:
        if os.path.exists(args.blast_results):
            blast_df = pd.read_csv(args.blast_results,sep="\t")
            blast_df.to_csv(f"{results_dir}/Full_Blast_results_against_PDB.tsv",sep="\t",index=False) #save it again in the new results folder
        else:
            blast_df = start_blast(args,results_dir)

    if args.pdb_results:
        if os.path.exists(args.pdb_results):
            print("File -Full_Blast_results_against_PDB_Filtered.tsv- found, reading and making a copy")
            PDB_files_dataframe = pd.read_csv(args.pdb_results,sep="\t")
            PDB_files_dataframe.to_csv(f"{results_dir}/Full_Blast_results_against_PDB_Filtered.tsv",sep="\t")
        else:
            PDB_files_dataframe = LASutils.process_blast_results(args, blast_df)
    else:
        PDB_files_dataframe = LASutils.process_blast_results(args,blast_df)


    #Read the file with paths to codeml output

    Codeml = pd.read_csv(os.path.join(args.codeml_output,"Gene_paths.txt"), sep='\s+', header=None)


    paml_results= Codeml.iloc[:,0].tolist()
    M8s_paths= Codeml.iloc[:,1].tolist() #paths
    dict_paths=dict(zip(paml_results,M8s_paths))

    fasta_dict = dict(zip(blast_df["Gene"],zip(blast_df["Prot_seq_clean"],blast_df["Prot_seq_full"])))


    for index,row in PDB_files_dataframe.iterrows():

        gene_name = row["Gene"]
        pdb_file_name =row["PDB_files"].strip("'")
        chains =','.join(list(row['Chains'].strip("'")))#I don't think I want a list here...

        prot_missing_data,prot_no_missing_data = fasta_dict[gene_name]
        prot_missing_data = prot_missing_data[0].replace("'","")  if isinstance(prot_missing_data,list) else prot_missing_data.replace("'","").replace("[","").replace("]","")
        prot_no_missing_data = prot_no_missing_data[0].replace("'","")  if isinstance(prot_no_missing_data,list) else prot_no_missing_data.replace("'","").replace("[","").replace("]","")
        filename_ent = os.path.join(args.pdb_files + "/pdb%s.ent" % pdb_file_name.lower())
        filename_cif = os.path.join(args.pdb_files + "/%s.cif" % pdb_file_name.lower())
        #if not os.path.exists(filename_ent) and not os.path.exists(filename_cif):
        if os.path.exists(filename_ent):
            residues_ID, pdb_sequence = LASutils.extract_sequence_from_PDB(filename_ent, chains,parser="ent")
        elif os.path.exists(filename_cif):
            residues_ID, pdb_sequence = LASutils.extract_sequence_from_PDB(filename_cif, chains,parser="cif")
        else:continue

        M8_file = ''.join(dict_paths[gene_name])
        if pdb_sequence:
            # Extract positions that possibly belong to a domain in index 1 format
            domains_idx = LASutils.PROSITE_domains(pdb_sequence)
        else:
            domains_idx = []

        sequences_dict = {"prot_missing_data":prot_missing_data,"prot_no_missing_data":prot_no_missing_data}

        if pdb_sequence:
            data = build_PyMOL_dataframe(pdb_sequence, sequences_dict, gene_name, chains, M8_file, domains_idx, residues_ID, pdb_file_name, args)


def visualize_single_result(args):
    """"""
    #single call to build dataframe and pymol
