import functools
from functools import partial
from Bio.Seq import Seq
from Bio.PDB import PDBList
from Bio import SeqRecord, SeqIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import MMCIFParser, MMCIF2Dict
import Bio.PDB as PDB
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1
from Bio.ExPASy import ScanProsite
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import pandas as pd
import os,shutil,ntpath,itertools
import itertools
import numpy as np


def folders(folder_name,basepath,overwrite=True):
    """ Creates a folder at the indicated location. It rewrites folders with the same name
    :param str folder_name: name of the folder
    :param str basepath: indicates the place where to create the folder
    """
    #basepath = os.getcwd()
    if not basepath:
        newpath = folder_name
    else:
        newpath = basepath + "/%s" % folder_name
    if not os.path.exists(newpath):
        try:
            original_umask = os.umask(0)
            os.makedirs(newpath, 0o777)
        finally:
            os.umask(original_umask)
    else:
        if overwrite:
            print("removing subdirectories") #if this is reached is because you are running the folders function twice with the same folder name
            shutil.rmtree(newpath)  # removes all the subdirectories!
            os.makedirs(newpath,0o777)
        else:
            pass
def url_reader(lines):
    """Matches a pattern and reads and stores the lines until it finds another identical pattern"""
    buffer = []
    for string in lines:
        if string.startswith('><a name'):
            if buffer: yield buffer
            buffer = [string]
        else:
            buffer.append(string)
    yield buffer
def fasta_to_sequences(File_Sequences,format):
    '''Same as before, can be merged with previous function'''

    # List_of_sequences = []
    # List_of_names = []
    fasta_sequences = SeqIO.parse(open(File_Sequences),format)
    # for fasta in fasta_sequences:
    #         name, sequence = fasta.id, str(fasta.seq)
    #         List_of_sequences.append(sequence)
    #         List_of_names.append(name)

    def extract_ids(fasta):
        return (fasta.id,str(fasta.seq))

    results = list(map(lambda fasta: extract_ids(fasta),fasta_sequences))

    zipped_results = list(zip(*results))

    id_names = zipped_results[0]
    sequences = zipped_results[1]


    return sequences,id_names
def validate(seq, alphabet='dna'):
    """
    Check that a sequence only contains values from DNA alphabet """
    import re
    alphabets = {'dna': re.compile('^[acgtn]*$', re.I),
             'protein': re.compile('^[acdefghiklmnpqrstvwy]*$', re.I)}
    if alphabets[alphabet].search(seq) is not None:
         return True
    else:
         return False
def translate_sequence(sequence,remove_missing_data = False): #Only for coding nucleotide sequences!!!!
    #from Bio.Alphabet import ProteinAlphabet
    if validate(sequence) == True: #Check if is a nucleotide sequence
        aa_seq = Seq(sequence).translate(stop_symbol="X")
    else:
        aa_seq=Seq(sequence)

    if remove_missing_data:
        aa_seq = Seq(str(aa_seq).replace('X',''))

    return aa_seq
def process_blast_results(args,blast_results):

    """Builds a dataframe with the following columns:

    [Gene]: Gene identifier,
    [PDB_files]: PDB files of homologous proteins
    [Chains]: Protein structure chain to which the protein is homologous
    [Percent ID]: Percent identity to the found homologous PDB structure
    [Resolution]: PDB file resolution, used for quality filtering
    [Prot_seq]: Translated sequence from the gene

    """

    print("Processing Blast results")



    # ungroup initial dataframe to have full control and do the next step
    s1 = blast_results.PDB_files.apply(str).str.strip('[]').str.replace(',', '').str.split(expand=True).apply(pd.Series, 1).stack()
    s2 = blast_results.Chains.apply(str).str.strip('[]').str.replace(',', '').str.split(expand=True).apply(pd.Series, 1).stack()
    s3 = blast_results.PercentID.apply(str).str.strip('[]').str.replace(',', '').str.split(expand=True).apply(pd.Series, 1).stack()
    s4 = blast_results.Prot_seq_full.apply(str).str.strip('[]').str.replace(',', '').replace("'","").str.split(expand=True).apply(pd.Series, 1).stack()
    #Resolving indexing
    s1.index=s1.index.droplevel(-1)
    s1.name='PDB_files'
    s2.index=s2.index.droplevel(-1)
    s2.name='Chains'
    s3.index=s3.index.droplevel(-1)
    s3.name='PercentID'
    s4.index=s4.index.droplevel(-1)
    s4.name='Prot_seq_full'
    # Join back the gene/protein columns
    ungrouped_dataframe = blast_results.drop(['PDB_files', 'Chains', 'PercentID','Prot_seq_full'], axis=1).join(s1)
    ungrouped_dataframe['Chains'] = s2
    ungrouped_dataframe['PercentID'] = s3
    ungrouped_dataframe['Prot_seq_full'] = s4

    ungrouped_dataframe.to_csv("{}/{}".format(args.results_dir,'Full_Blast_results_against_PDB.tsv'), sep='\t', index=False)
    # change the dtype of %ID and coverage columns
    ungrouped_dataframe[['PercentID']] = ungrouped_dataframe[['PercentID']].apply(pd.to_numeric)
    #Reduce the amount of PDB files: First, higher Identity. Get 3 first values
    ungrouped_dataframe=ungrouped_dataframe.groupby(['Gene'],as_index=False)[['Gene','PDB_files','Chains','PercentID','Prot_seq_full']].apply(lambda x: x.nlargest(args.number_homologous, columns=['PercentID']))
    ungrouped_dataframe.reset_index(inplace=True,drop=True)

    #Download the selected PDB files :
    pdbl = PDBList()
    parser = PDBParser(PERMISSIVE=1)
    ungrouped_dataframe['Resolution'] = '' #initialize an empty column

    for index, row in ungrouped_dataframe.iterrows():
        pdb_file_name = row["PDB_files"].strip("'")
        filename_ent = os.path.join(args.pdb_files + "/pdb%s.ent" % pdb_file_name.lower())
        filename_cif = os.path.join(args.pdb_files + "/%s.cif" % pdb_file_name.lower())
        if not os.path.exists(filename_ent) and not os.path.exists(filename_cif):
            pdbl.retrieve_pdb_file(pdb_file_name, pdir=args.pdb_files)
        else:
            print("PDB file found, reading")
            if os.path.exists(filename_ent):
                structure = parser.get_structure(row["PDB_files"], filename_ent)
                resolution = float(structure.header["resolution"])
            else: #.cif files
                mmcif_dict = MMCIF2Dict.MMCIF2Dict(filename_cif)
                if "_reflns.d_resolution_high" in mmcif_dict.keys():
                    resolution = float(mmcif_dict["_reflns.d_resolution_high"][0])
                else:
                    resolution= float(5)
            ungrouped_dataframe.loc[index, 'Resolution'] = resolution

    # def retrieve_pdb(row,ungrouped_dataframe):
    #     pdb_file_name = row["PDB_files"].strip("'")
    #     filename_ent = os.path.join(args.pdb_files + "/pdb%s.ent" % pdb_file_name.lower())
    #     filename_cif = os.path.join(args.pdb_files + "/%s.cif" % pdb_file_name.lower())
    #     if not os.path.exists(filename_ent) and not os.path.exists(filename_cif):
    #         pdbl.retrieve_pdb_file(pdb_file_name, pdir=args.pdb_files)
    #     else:
    #         print("PDB file found, reading")
    #         if os.path.exists(filename_ent):
    #             structure = parser.get_structure(row["PDB_files"], filename_ent)
    #             resolution = float(structure.header["resolution"])
    #         else: #.cif files
    #             mmcif_dict = MMCIF2Dict.MMCIF2Dict(filename_cif)
    #             if "_reflns.d_resolution_high" in mmcif_dict.keys():
    #                 resolution = float(mmcif_dict["_reflns.d_resolution_high"][0])
    #             else:
    #                 resolution= float(5)
    #         ungrouped_dataframe.loc[index, 'Resolution'] = resolution

    #results = list(map(lambda row: partial(retrieve_pdb(row,ungrouped_dataframe)), ungrouped_dataframe.iterrows())) #TODO: Think about


    ungrouped_dataframe['Resolution'] = ungrouped_dataframe['Resolution'].astype('float64')
    filtered_dataframe = ungrouped_dataframe.groupby(['Gene'],as_index=False)[['Gene','PDB_files','Chains','PercentID','Resolution','Prot_seq_full']].apply(lambda x: x.nsmallest(int(args.number_homologous), columns=['Resolution']))
    filtered_dataframe.to_csv('{}/Full_Blast_results_against_PDB_Filtered.tsv'.format(args.results_dir),sep='\t',index=False)

    return filtered_dataframe
def extract_sequence_from_PDB(PDB_file,chains,parser="ent"):
    ''' Returns both the sequence contained in the PDB file and the residues number identifiers for the desired chains'''
    name = ntpath.basename(PDB_file).split('.')[0]
    if parser == "ent":
        parser = PDB.PDBParser()
        structure = parser.get_structure('%s' % (name), PDB_file)
    elif parser == "cif":
        parser = MMCIFParser()
        structure = parser.get_structure('%s' % (name), PDB_file)

    #chains='all'
    ############## Iterating over residues to extract all of them even if there is more than 1 chain
    sequence = []
    residues_ID = []
    chains_found = 0
    if chains == 'all':
        for chain in structure.get_chains():
            for residue in chain:
                if is_aa(residue.get_resname(), standard=True):
                    sequence.append(residue.get_resname())
                    residues_ID.append(residue.get_id()[1])
    else :
        accumulated_residues = []
        accumulated_ids = []
        for chain_letter in chains:
           if chain_letter in structure[0].child_dict.keys():
                for residue in structure[0][chain_letter]:
                    if is_aa(residue.get_resname(), standard=True):
                        accumulated_residues.append(residue.get_resname())
                        accumulated_ids.append(residue.get_id()[1])
                chains_found += 1
        sequence.append(''.join(accumulated_residues))
        residues_ID.append(accumulated_ids)

    joined_sequence = ''.join(sequence)
    PDB_sequence = seq1(joined_sequence) #3 letter code into 1 letter code

    if chains_found >= 1:
        residues_ID = list(itertools.chain.from_iterable(residues_ID)) #join list of lists


    return residues_ID,PDB_sequence
def PROSITE_domains(PDB_sequence):
    "Returns domain positions found in the PDB sequence in 1-index format"
    handle = ScanProsite.scan(seq=PDB_sequence,mirror="https://prosite.expasy.org/") #can use a sequence (seq=) or a pdbID, but we want a specific chain
    result = ScanProsite.read(handle)
    domains=[]
    for index,segment in enumerate(result,1):
        #if 'score' in segment:
            domains.append(list(range(segment['start'],segment['stop'])))
        #else:
            #pass

    if len(domains) >= 1:
        return list(itertools.chain.from_iterable(domains))
    else:
        return domains
def extract_positive_sites(file, prob):  # file = Out file from M8 #prob = 95 or 99
    ''' Reads M8 from M7vsM8 Codeml test: Will return the Positions of the selected sites with regards to the ALIGNMENT, later use function to find the equivalent positions in the clean from missing data sequence'''
    length_alignment = []
    #/home/dragon/drive/lys/Dropbox/MasterRepositories/LYS_Automatic_Search/TestSequences/Fasta_sca1.106.1.map2per.fasta.paml/M8/out
    with open(file, 'r') as f:
        data = f.read().split('\n')
        positions = []
        line_start = []
        line_stop = []  # Line number reference to stop collecting the info related to positive selected sites
        for number, line in enumerate(data, 1):
            if 'Bayes Empirical Bayes (BEB) analysis (Yang, Wong & Nielsen 2005. Mol. Biol. Evol. 22:1107-1118)' in line:
                line_start = number
            if 'The grid' in line:
                line_stop = number
        if 'Bayes Empirical Bayes (BEB) analysis (Yang, Wong & Nielsen 2005. Mol. Biol. Evol. 22:1107-1118)' in data:
            diff = line_stop - line_start - 3
            for i in range(6, int(diff)):  # Start at line 6 after match
                position = data[data.index(
                    'Bayes Empirical Bayes (BEB) analysis (Yang, Wong & Nielsen 2005. Mol. Biol. Evol. 22:1107-1118)') + i]
                if prob == 99:
                    if str(position.split()[2]).endswith('**'):  # * > 95% confidence ** > 99 confidence
                        # print(position)
                        position = position.split()[0]
                        positions.append(int(position))  # Needs to be an int for later comparison
                else:
                    if str(position.split()[2]).endswith('*'):  # * > 95% confidence ** > 99 confidence
                        # print(position)
                        position = position.split()[0]
                        positions.append(int(position))  # Needs to be an int for later comparison
        return positions
def global_alignment(chain_A,chain_B):
    '''Global alignment between 2 given chains(str format)'''
    print("Performing global alignment")
    alignment_global = pairwise2.align.globalms(chain_A, chain_B, 2, -1, -5, -.1)
    #print(format_alignment(*alignment_global[0])) #print the alignment if desired
    print(alignment_global)

    alignment_info_global = alignment_global[0]
    #Aligned_3u84, Aligned_ENST00000337652, score, begin, end = alignment_info_global
    return alignment_info_global,alignment_global


def equivalent_positions(chain_a, chain_b,aligned_a, aligned_b,residues_ID = None, domain_positions =None): #chainA = PDB #chainB=Transcript #Positions = non_gapped_sites of Positive Selected sites (optional argument because we will use this function in 2 manners)
    ''' This function returns the corresponding coordinates of the gene residues in the PDB sequence, the domain positions in the PDB sequence
    or among the gene with missing data and without
    :param chain_a
    :param chain_b
    :param aligned_a
    :param aligned_b

    '''
    print("Mapping positions")

    #TODO: Probably it can be better managed with a dictionary?
    #OriginalA = np.array(range(1, len(chain_a))) #Array of positions of the first given sequence

    # Find the non-gapped equivalent positions
    non_gapped_sites = [] #Will store the Positions where there are no gaps in the aligned sequences, these are still not the PDB positions !!
    # for index,residue in enumerate(aligned_a):
    #     if aligned_a[index] != '-' and aligned_b[index] != '-':
    #             #OriginalA[:index] += 1
    #             non_gapped_sites.append(index +1) #In index 1 for PDB file
    #     else:
    #         pass

    non_gapped_sites = list(map(functools.partial(lambda aligned_a,aligned_b,index: index +1 if aligned_a[index] != '-' and aligned_b[index] != '-' else None, aligned_a,aligned_b),list(range(len(aligned_a)))))
    non_gapped_sites=list(filter(lambda v: v is not None, non_gapped_sites))

    final_positions = ['nan']*len(chain_b) #should have the lenght of the protein. If equivalence is found nan will be replaced with the equivalent PDB positions
    #Highlight: Finding the equivalent PDB Position of the residues in the gene transcript ignoring the GAPS generated in the alignment#######
    positions_in_PDB =[]
    gaps_first_segment = ''.join(aligned_a[0:non_gapped_sites[0]]).count('-')
    positions_in_PDB.append(non_gapped_sites[0] - gaps_first_segment)

    #Highlight: Finding out the equivalent positions in the gene transcript of the PDB residues######: Same thing the other way round
    positions_in_transcript = []
    gaps_first_segment_transcript = ''.join(aligned_b[0:non_gapped_sites[0]]).count('-')
    positions_in_transcript.append(non_gapped_sites[0] - gaps_first_segment_transcript)


    accumulated_number_gaps_in_this_segment = gaps_first_segment
    accumulated_number_gaps_in_this_segment_transcript = gaps_first_segment_transcript
    for i in range(0,len(non_gapped_sites)-1): #cannot start at position 1 otherwise we skip one
        #print(i)
        #try:
        accumulated_number_gaps_in_this_segment += ''.join(aligned_a[non_gapped_sites[i]:non_gapped_sites[i+1]]).count('-')
        positions_in_PDB.append(non_gapped_sites[i+1] - accumulated_number_gaps_in_this_segment)
        # except:
        #     pass
        #try:
        accumulated_number_gaps_in_this_segment_transcript += ''.join(aligned_b[non_gapped_sites[i]:non_gapped_sites[i + 1]]).count('-')
        positions_in_transcript.append(non_gapped_sites[i+1] - accumulated_number_gaps_in_this_segment_transcript) # plus on otherwise negative numbers
        # except:
        #     pass

    #Position_in_Transcript.append(non_gapped_sites[1] - gaps_first_segment_Transcript - ''.join(aligned_b[non_gapped_sites[0]:non_gapped_sites[1]]).count('-'))

    #
    # for i in range(0, len(non_gapped_sites)):  # we skip the first interval
    #     try:
    #         accumulated_number_gaps_in_this_segment_transcript += ''.join(aligned_b[non_gapped_sites[i]:non_gapped_sites[i + 1]]).count('-')
    #         position_in_transcript.append(non_gapped_sites[i+1] - accumulated_number_gaps_in_this_segment_transcript) # plus on otherwise negative numbers
    #     except:
    #         pass
    if not domain_positions:
        if residues_ID: #residues ID has the same length as the PDB crystallized sequence #recall: is cheaper to do the if statement outside once
            for position_transcript,position_PDB in zip(positions_in_transcript,positions_in_PDB): #this are non_gapped_sites of the correspondance on both sides
                    final_positions[position_transcript - 1] = residues_ID[position_PDB - 1]
        else:
            for position_transcript, position_PDB in zip(positions_in_transcript,positions_in_PDB):  # this are non_gapped_sitess of the correspondance on both sides
                    final_positions[position_transcript - 1] = position_PDB

        return final_positions
    else:
        for position_transcript, position_PDB in zip(positions_in_transcript, positions_in_PDB):
            final_positions[position_transcript - 1] = position_PDB
        equivalent_domain_positions = []
        for residue_id, domain_position in zip( residues_ID,domain_positions): #Finding the correspondant coordinates of the functional domains in the PDB file structure
            try: #TODO: Think How to remove try/except statemts
                specific_corresponding_domain_position = final_positions[domain_position - 1]
                try:
                    equivalent_domain_positions.append(residues_ID[specific_corresponding_domain_position - 1])
                except:
                    pass
            except:
                pass
        return equivalent_domain_positions

def corresponding_positions_missing_notmissing_data(prot_missing_data,prot_no_missing_data): #Missing_data = WITH missing data ; Clean = NO missing data
    ''' Returns list of the equivalent positions among 2 sequences, in this case the same sequence with and without missing data'''

    alignment_info_global, alignment_global = global_alignment(prot_missing_data,prot_no_missing_data)
    aligned_A, aligned_B, score, begin, end = alignment_info_global #Sequence A = with missing data ('N'); Sequence B = Not missing data
    mapped_positions = equivalent_positions(prot_missing_data,prot_no_missing_data, aligned_A, aligned_B) #Map corresponding positions of the sequence with missing data to the sequence without missing data
    print("---------------------------------")
    return mapped_positions

def corresponding_coordinates_and_labels_PDB_gene(pdb_sequence,pdb_file_name,prot_gene,positive_sites,domains_idx,residues_ID,basepath,args,gene_name,mapped_positions = None):
    #pdb_sequence,prot_no_missing_data, positives_sites,domains_idx,residues_ID,args.proteins,args.print_alignment,gene_name, mapped_positions

    '''Performs the alignments, retrieves the equivalent coordinates and labels among the PDB sequence and the gene for each of the PDB positions.
    Produces a dataframe with different labels for each position, where 'Not'refers to not being positively selected and not functional domain,
    'Selected' stands for positive selected, 'Domain', states that it belongs to the functional domain or 'Selected_and_Domain' which stands for both
    :param: list domains_idx: Prosites results
    '''

    prot_gene = str(prot_gene)
    #Perform the global alignment
    alignment_info_global, alignment_global = global_alignment(pdb_sequence, prot_gene)
    aligned_a, aligned_b, score, begin, end = alignment_info_global #A = PDB; B = Gene(with or without missing data)
    if args.print_alignment == 'yes':
        print(format_alignment(*alignment_global[0]))

    prot_gene_positions = list(range(1, len(prot_gene) + 1)) #create a list of numbers in range 1 to gene length
    #Extract the corresponding positions of the positive selected sites in the clean of missing data gene sequence
    positive_positions = []
    if mapped_positions: #If the gene sequence has been cleaned (removed missing information) and we have the corresponding mapped positions
        for element in positive_sites:
            try:
                positive_positions.append(mapped_positions.index(element)) #Append the index of the residue that has the same id as the one in positive selected sites
            except:
                pass
    else: #we don't do anything
        positive_positions = positive_sites


    #For the dataframe we can label the positions from 1 to length of gene sequence (with or without missing data)
    positions_dataframe = pd.DataFrame(pd.Series(prot_gene_positions))
    positions_dataframe.rename(columns={positions_dataframe.columns[0]: "Gene_Position"}, inplace=True)


    mapped_positions = equivalent_positions(pdb_sequence, prot_gene, aligned_a, aligned_b,residues_ID = residues_ID) #List of equivalent positions of each of the gene residues into the PDB structure(residues ID!)
    ###The functional domains require also processing to find their equivalent residues in the PDB
    domains_idx = [residues_ID[index-1] for index in domains_idx ] #Index 1 in both, residues ID and Prosite,so fix for python 0 indexing

    #Add the positions to the dataframe
    positions_dataframe['PDB_Position'] = mapped_positions
    Label_1 = positions_dataframe['Gene_Position'].isin(positive_positions) #Check if the position in the gene sequence is + selected
    Label_2 = positions_dataframe['PDB_Position'].isin(domains_idx) #Check if the PDB position is a functional domain

    positions_dataframe['Label_1'] = Label_1  # Create a Column were the positive positions have a label True
    positions_dataframe['Label_2'] = Label_2 #Create a column were the functional domain positions recieve the label True
    positions_dataframe['Label_1'] = positions_dataframe['Label_1'].replace(True, 'Selected')
    positions_dataframe['Label_1'] = positions_dataframe['Label_1'].replace(False, 'Not')
    positions_dataframe['Label_2'] = positions_dataframe['Label_2'].replace(True, 'Domain')
    positions_dataframe['Label_2'] = positions_dataframe['Label_2'].replace(False, 'Not')
    positions_dataframe['Label'] = pd.Series(prot_gene_positions)
    for index, row in positions_dataframe.iterrows():

        if positions_dataframe.iloc[index, positions_dataframe.columns.get_loc('Label_1')] == 'Selected' and positions_dataframe.iloc[index, positions_dataframe.columns.get_loc('Label_2')] == 'Domain':
            positions_dataframe.iloc[index, positions_dataframe.columns.get_loc('Label')] = 'Selected_and_Domain'
        elif positions_dataframe.iloc[index, positions_dataframe.columns.get_loc('Label_1')] == 'Selected':
            positions_dataframe.iloc[index, positions_dataframe.columns.get_loc('Label')] = 'Selected'
        elif positions_dataframe.iloc[index, positions_dataframe.columns.get_loc('Label_2')] == 'Domain':
            positions_dataframe.iloc[index, positions_dataframe.columns.get_loc('Label')] = 'Domain'
        else:
            positions_dataframe.iloc[index, positions_dataframe.columns.get_loc('Label')] = 'Not'
    positions_dataframe.drop(columns=['Label_1', 'Label_2'], axis=1, inplace=True)

    positions_dataframe.to_csv(args.results_dir + "/%s_%s_Positions.tsv" % (gene_name,pdb_file_name), index=False, sep='\t')
    print('Dataframe Ready at %s!' % (args.results_dir))
    return positions_dataframe





