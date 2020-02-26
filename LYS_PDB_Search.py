 #!/usr/bin/env python
#####IMPORT ALL NECESSARY MODULES#####

#CITATION https://www.biorxiv.org/content/10.1101/540229v1

import sys, ast, json
import os
import re
from urllib import request
import argparse
import itertools
import ntpath
import numbers
import decimal
import sys, os
#import pymol
#from pymol import cmd
#from pymol.cgo import *
#from pymol.vfont import plain
import ntpath
import pandas as pd #Windows users: I copy-paste the pandas,dateutil and pyltz folders from anaconda 2!! into the site-packages folder of pymol(only for temporary use, other wise it gets confused with the paths of the packages)
import numpy as np
from pandas import Series
from math import isnan
#Biopython
from Bio.Alphabet.IUPAC import ExtendedIUPACProtein
from Bio.Alphabet import generic_protein
from Bio import SeqRecord,Alphabet,SeqIO
from Bio.SeqIO import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import Bio.PDB as PDB
from Bio.Seq import MutableSeq
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.pairwise2 import format_alignment
from Bio import BiopythonWarning
import warnings
warnings.simplefilter('ignore', BiopythonWarning)
#Readline
try:
  import readline #Linux
except ImportError:
  import pyreadline as readline #Windows
import rlcompleter


#Command Line Arguments
parser = argparse.ArgumentParser()

#COMPULSORY files: Give error if not present
parser.add_argument('--Proteins',required = True,help = 'FULL Path to file with all positively selected peptide sequences')
parser.add_argument('--Codeml',required =True,help ='FULL path to tab separated text file with the fasta ids and the paths to codeml output files (M8)')
#OPTIONAL arguments:
parser.add_argument('--format', help = 'Sequence or Multiple Alignment File Format (default fasta), do not write it with quotes',default = 'fasta')
parser.add_argument("--prob", help = 'Choice of level of posterior probability on the sites, 95% or 99% from M8 Out file', default= int(99)) 	# prob : Posterior probability percentage of the positive sites , 99 or 95, default= 99
parser.add_argument("--missing_data", help = 'Decide if the missing data("N") should be removed from the nucleotide sequence, default = yes, recommended for the alignment', default='yes') 	# missing_data : remove or not missing data from the sequence alignment: yes or no, default=yes
parser.add_argument("--print_alignment",help= 'Choose to visualize the PDB file sequence aligned with the gene, default = no', default='no')

args = parser.parse_args()
Genes,codeml_output,Gene_file_format,prob,missing_data,print_alignment = [args.Proteins,args.Codeml,args.format,args.prob,args.missing_data,args.print_alignment]

#Check all the necessary files are present and the right arguments are given:
File_formats = ['fasta','phylip_sequential','clustal','embl','genebank','gb','abi','ace','fastq-sanger','fastq','fastq-solexa','fastq-illumina','ig','imgt','pdb-seqres','pdb-atom','phd','phylip','pir','seqxml','sff','stockholm','swiss','tab','qual','uniprot-xml']
if not Genes or not codeml_output:
    parser.error('Missing one of the 2 required arguments: Genes sequences / Dataframewith codeml paths and gene name')
    #Genes = "/home/lys/Downloads/LYS_Automatic_Search-master/TestSequences/Test5.fasta"
    #codeml_output = "/home/lys/Downloads/LYS_Automatic_Search-master/TestSequences/Fasta_paths.txt"
elif prob not in [95,99]:
    parser.error('Probability values can only be 95  or 99')
elif missing_data not in ['no','yes']:
    parser.error('Choose to keep missing data: yes or to remove it: no')
elif Gene_file_format not in File_formats:
    parser.error('Invalid Gene File Format.Check available SeqIO biopython file formats ')
else:
    print('Building dataframe...')

basepath = os.path.dirname(Genes)

def URL_reader(lines):
    """Matches a pattern and reads and stores the lines until it finds another identical pattern"""
    buffer = []
    for string in lines:
        if string.startswith('><a name'):
            if buffer: yield buffer
            buffer = [string]
        else:
            buffer.append(string)
    yield buffer
def Blast_and_List_PDB_files(gene_sequence): #Translated gene sequence
    """Using the protein sequences, blast them against the PDB database and parse the result, only significant e-values: 1e-6"""
    url = "http://www.rcsb.org/pdb/rest/getBlastPDB1?sequence=%s&eCutOff=1e-6&matrix=BLOSUM62" % (gene_sequence) #Reduce cut-off to 0
    PDB_IDs=[]
    Chains=[]
    Identities=[]
    Positives=[]
    data = request.urlopen(url).read().splitlines(True)
    data =[line.decode("utf-8") for line in data]
    for index,heading_and_lines in enumerate(URL_reader(data)):
        if index == 0:
            pass
        else:
            heading = heading_and_lines[0] #line with ><a name
            lines = heading_and_lines[4] #line with identities and positives
            try:
                if 'Identities' in lines:#avoids weird lines that can mess up the information extraction
                    #Headings
                    heading = re.sub(r'(><a name = [0-9]+></a>)', '', heading).strip() #Removes the first part of the string
                    PDB_IDs.append(heading.split(':')[0])
                    Chains.append(heading.split(':')[2].split('|')[0])
                    #Other lines
                    identity=lines.split(',')[0]
                    identity=eval(identity[identity.find("(")+1:identity.find(")")].strip('%'))
                    Identities.append(identity)
                    positives = lines.split(',')[1]
                    positives = eval(positives[positives.find("(") + 1:positives.find(")")].strip('%'))
                    Positives.append(positives)
                else:
                    pass
            except:
                pass
    return PDB_IDs,Chains, Identities,Positives
def Initial_dataframe(List_protein_sequences,List_names):
    """Builds a dataframe with all the pertinent information from the blast results"""
    dictionary ={}
    for gene,name in zip(List_protein_sequences,List_names):
        PDBs_List,Chains,PercentID,Coverage = Blast_and_List_PDB_files(gene)
        if not PDBs_List:#to avoid empty information
            pass
        else:
            dictionary[name] = PDBs_List,Chains,PercentID,Coverage
    df = pd.DataFrame.from_dict(dictionary,orient='index')
    df.index.name = 'Genes'
    df.reset_index(inplace=True)
    df.columns = ['Gene','PDB_files','Chains','PercentID','Coverage']
    #df.to_csv(os.path.join(basepath + 'Full_Blast_results_against_PDB.tsv'),sep='\t',index=False)
    return df
def Extract_sequence_from_PDB(PDB_file,chains):
    ''' Returns both the sequence contained in the PDB file and the residues coordinates for the desired chains'''
    from Bio.PDB.PDBParser import PDBParser
    from Bio.PDB import MMCIFParser
    Name = ntpath.basename(PDB_file).split('.')[0]
    try:
        parser = PDB.PDBParser()
        structure = parser.get_structure('%s' % (Name), PDB_file)
    except:
        parser = MMCIFParser()
        structure = parser.get_structure('%s' % (Name), PDB_file)

    #chains='all'
    ############## Iterating over residues to extract all of them even if there is more than 1 chain
    sequence = []
    Residues_ID = []
    if chains == 'all':
        for chain in structure.get_chains():
            for residue in chain:
                if is_aa(residue.get_resname(), standard=True):
                    sequence.append(residue.get_resname())
                    Residues_ID.append(residue.get_id()[1])
    else :
        accumulated_residues = []
        accumulated_ids = []
        for letter in chains:
           try: #in case the chain requested does not exits
                for residue in structure[0][letter]:
                    #print(letter)
                    if is_aa(residue.get_resname(), standard=True):
                        accumulated_residues.append(residue.get_resname())
                        accumulated_ids.append(residue.get_id()[1])
           except:
               pass
        sequence.append(''.join(accumulated_residues))
        Residues_ID.append(accumulated_ids)

    joined_sequence = ''.join(sequence)
    PDB_sequence = seq1(joined_sequence) #3 letter code into 1 letter code
    try:
        Residues_ID = list(itertools.chain.from_iterable(Residues_ID))
    except:
        pass

    return Residues_ID,PDB_sequence
def Download_selected_PDB_files_and_add_resolution_and_percentage_id(initial_dataframe):
    import os
    initial_dataframe.to_csv(os.path.join(basepath + 'Initial_dataframe.tsv'),sep='\t',index=False)
    # ungroup initial dataframe to have full control and do the next step
    s1 = initial_dataframe.PDB_files.apply(str).str.strip('[]').str.replace(',', '').str.split(expand=True).apply(Series, 1).stack()
    s2 = initial_dataframe.Chains.apply(str).str.strip('[]').str.replace(',', '').str.split(expand=True).apply(Series, 1).stack()
    s3 = initial_dataframe.PercentID.apply(str).str.strip('[]').str.replace(',', '').str.split(expand=True).apply(Series, 1).stack()
    s4 = initial_dataframe.Coverage.apply(str).str.strip('[]').str.replace(',', '').str.split(expand=True).apply(Series, 1).stack()
    #Resolving indexing
    s1.index=s1.index.droplevel(-1)
    s1.name='PDB_files'
    s2.index=s2.index.droplevel(-1)
    s2.name='Chains'
    s3.index=s3.index.droplevel(-1)
    s3.name='PercentID'
    s4.index=s4.index.droplevel(-1)
    s4.name='Coverage'

    #Join back the gene columns

    ungrouped_dataframe = initial_dataframe.drop(['PDB_files', 'Chains', 'PercentID', 'Coverage'], axis=1).join(s1)
    ungrouped_dataframe['Chains'] = s2
    ungrouped_dataframe['PercentID'] = s3
    ungrouped_dataframe['Coverage'] =s4

    ungrouped_dataframe.to_csv(os.path.join(basepath + 'Full_Blast_results_against_PDB.tsv'),sep='\t',index=False)
    #change the dtype of %ID and coverage columns
    ungrouped_dataframe[['PercentID', 'Coverage']] = ungrouped_dataframe[['PercentID', 'Coverage']].apply(pd.to_numeric)
    #ungrouped_dataframe = initial_dataframe
    #Reduce the amount of PDB files: First, higher Identity, second, coverage. Get 3 first values
    max_vals=ungrouped_dataframe.groupby(['Gene'])[['PDB_files','Chains','PercentID','Coverage']].apply(lambda x: x.nlargest(3, columns=['PercentID']))
    indexes =max_vals.index.levels[1]
    max_vals=ungrouped_dataframe.groupby(['Gene'])[['PDB_files','Chains','PercentID','Coverage']].apply(lambda x: x.nlargest(3, columns=['Coverage']))
    indexes = max_vals.index.levels[1]
    ungrouped_dataframe=ungrouped_dataframe.iloc[indexes] #Get only the rows where the percentid and coverage are maximum
    ungrouped_dataframe.reset_index(inplace=True,drop=True)

    #Donwload the selected PDB files :
    from Bio.PDB import PDBList
    import os
    pdbl = PDBList()
    from Bio.PDB.PDBParser import PDBParser
    from Bio.PDB import MMCIFParser, MMCIF2Dict
    parser = PDBParser(PERMISSIVE=1)
    ungrouped_dataframe['Resolution'] = ''
    for index, row in ungrouped_dataframe.iterrows():
        try:
            if not basepath:
                pdbl.retrieve_pdb_file(row["PDB_files"].strip("'"), pdir='PDB_files')
            else:
                pdbl.retrieve_pdb_file(row["PDB_files"].strip("'"),pdir= os.path.join(basepath + '/PDB_files'))

        except:
            pass
        #Resolution
        if not basepath: #If path is in the same directory
            try: #.ent files
                filename = os.path.join('PDB_files'+ "/pdb%s.ent" % row["PDB_files"].strip("'").lower())
                structure = parser.get_structure(row["PDB_files"], filename)
                resolution = structure.header["resolution"]
            except: #.cif files
                filename = os.path.join('PDB_files'+ "/%s.cif" % row["PDB_files"].strip("'").lower())
                mmcif_dict = MMCIF2Dict.MMCIF2Dict(filename)
                try:
                    resolution = mmcif_dict["_reflns.d_resolution_high"]
                except:
                    resolution= np.nan


        else:
            try:
                filename = os.path.join(basepath + '/PDB_files' + "/pdb%s.ent" % row["PDB_files"].strip("'").lower())
                structure = parser.get_structure(row["PDB_files"], filename)
                resolution = structure.header["resolution"]
            except:
                filename = os.path.join(basepath + '/PDB_files' + "/%s.cif" % row["PDB_files"].strip("'").lower())
                mmcif_dict = MMCIF2Dict.MMCIF2Dict(filename)
                try:
                    resolution = mmcif_dict["_reflns.d_resolution_high"]
                except:
                    resolution = np.nan

        #structure = parser.get_structure(row["PDB_files"], filename)
        #resolution = structure.header["resolution"]
        try:
            ungrouped_dataframe.loc[index,'Resolution'] = float(resolution)
        except:#assign value to missing resolution values to force to discard them
            ungrouped_dataframe.loc[index, 'Resolution']= float(5)
    ungrouped_dataframe['Chains'] = ungrouped_dataframe['Chains']

    ungrouped_dataframe['Resolution'] = ungrouped_dataframe['Resolution'].astype('float64')

    min_vals = ungrouped_dataframe.groupby(['Gene'])[['PDB_files','Chains','PercentID','Coverage','Resolution']].apply(lambda x: x.nsmallest(3, columns=['Resolution']))
    min_indexes = min_vals.index.levels[1]
    filtered_dataframe = ungrouped_dataframe.iloc[min_indexes]  # Get only the rows where the percentid and coverage are maximum
    filtered_dataframe.reset_index(inplace=True, drop=True)
    filtered_dataframe.to_csv('Full_Blast_results_against_PDB_Filtered.tsv',sep='\t',index=False)
    return filtered_dataframe
def Global_alignment(chain_A,chain_B):
    '''Global alignment between 2 given chains(str format)'''
    from Bio import pairwise2
    from Bio.pairwise2 import format_alignment
    alignment_global = pairwise2.align.globalms(chain_A, chain_B, 2, -1, -5, -.1)
    #print(format_alignment(*alignment_global[0])) #print the alignment if desired
    alignment_info_global = alignment_global[0]
    Aligned_3u84, Aligned_ENST00000337652, score, begin, end = alignment_info_global
    return alignment_info_global,alignment_global
def Local_alignment(chain_A,chain_B):
    '''Local alignment between 2 given chains'''
    from Bio import pairwise2
    from Bio.pairwise2 import format_alignment
    alignment_local = pairwise2.align.localms(chain_A, chain_B, 2, -1, -5, -.1) ##0.5 points are deducted when opening a gap, and 0.1 points are deducted when extending it.
    alignment_info_local = alignment_local[0]
    Aligned_A, Aligned_B, score, begin, end = alignment_info_local
    return alignment_info_local,alignment_local
def fasta_to_sequences(File_Sequences,Format):
    '''Same as before, can be merged with previous function'''
    from Bio import SeqRecord, SeqIO
    List_of_sequences = []
    List_of_names = []
    fasta_sequences = SeqIO.parse(open(File_Sequences),Format)
    for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            List_of_sequences.append(sequence)
            List_of_names.append(name)
    return List_of_sequences,List_of_names
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
def Translate_sequence(sequence): #Only for coding nucleotide sequences!!!!
    from Bio.Alphabet import IUPAC,ProteinAlphabet

    if validate(sequence) == True: #Check if is a nucleotided sequence
        aa_seq = Seq(sequence).translate(stop_symbol="X")
    else:
        aa_seq=Seq(sequence)

    return aa_seq
def Translate_and_Remove_missing_data(sequence):
    '''Remove the missing data ('X') from the sequences after being translated, otherwise the codons are affected'''

    if validate(sequence) == True:
        clean_sequence = Seq(sequence).translate(stop_symbol="X")
        clean_sequence = clean_sequence.ungap('X')
    else:
        clean_sequence = Seq(sequence).ungap('X') #maybe ungap

    return clean_sequence
def equivalent_positions(chain_A, chain_B,Aligned_A, Aligned_B,Residues_ID = None, Domain_Positions =None): #chainA = PDB #chainB=Transcript #Positions = List of Positive Selected sites (optional argument because we will use this function in 2 manners)
    ''' This function returns the corresponding coordinates of the gene residues in the PDB sequence, the domain positions in the PDB sequence
    or among the gene with missing data and without'''
    import numpy as np
    #OriginalA = np.array(range(1, len(chain_A))) #Array of positions of the first given sequence

    ##Equivalent positions
    List = [] #Will store the Positions where there are no gaps in the aligned sequences, these are still not the PDB positions !!
    for index,residue in enumerate(Aligned_A):
        if Aligned_A[index] != '-' and Aligned_B[index] != '-':
                #OriginalA[:index] += 1
                List.append(index +1) #In index 1 for PDB file
        else:
            pass
    Final_Positions = ['nan']*len(chain_B) #should have the lenght of the gene. If equivalence is found nan will be replaced with the equivalent PDB positions
    #####Finding the equivalent PDB Position of the residues in the gene transcript ignoring the GAPS generated in the alignment#######
    Position_in_PDB =[]
    gaps_first_segment = ''.join(Aligned_A[0:List[0]]).count('-')
    Position_in_PDB.append(List[0] - gaps_first_segment)
    #Position_in_PDB.append(List[1] - gaps_first_segment - ''.join(Aligned_A[List[0]:List[1]]).count('-'))

    accumulated_number_gaps_in_this_segment = gaps_first_segment
    for i in range(0,len(List)): #we skip the first interval, already performed !
        try:
            accumulated_number_gaps_in_this_segment += ''.join(Aligned_A[List[i]:List[i+1]]).count('-')
            Position_in_PDB.append(List[i+1] - accumulated_number_gaps_in_this_segment)
        except:
            pass
    #####Finding out the equivalent positions in the gene transcript of the PDB residues######: Same thing the other way round
    Position_in_Transcript = []
    gaps_first_segment_Transcript = ''.join(Aligned_B[0:List[0]]).count('-')
    Position_in_Transcript.append(List[0] - gaps_first_segment_Transcript)
    #Position_in_Transcript.append(List[1] - gaps_first_segment_Transcript - ''.join(Aligned_B[List[0]:List[1]]).count('-'))

    accumulated_number_gaps_in_this_segment_transcript = gaps_first_segment_Transcript
    for i in range(0, len(List)):  # we skip the first interval
        try:
            accumulated_number_gaps_in_this_segment_transcript += ''.join(Aligned_B[List[i]:List[i + 1]]).count('-')
            Position_in_Transcript.append(List[i+1] - accumulated_number_gaps_in_this_segment_transcript) # plus on otherwise negative numbers
        except:
            pass
    Equivalent_Domain_positions = []
    if not Domain_Positions:
        for position_transcript,position_PDB in zip(Position_in_Transcript,Position_in_PDB): #this are lists of the correspondance on both sides

                 if Residues_ID: #residues ID has the same length as the PDB crystallized sequence
                    Final_Positions[position_transcript - 1] = Residues_ID[position_PDB - 1]
                 else:
                    Final_Positions[position_transcript-1]  = position_PDB

        return Final_Positions
    else:
        for position_transcript, position_PDB in zip(Position_in_Transcript, Position_in_PDB):
            Final_Positions[position_transcript - 1] = position_PDB

        for residue_id, domain_position in zip( Residues_ID,Domain_Positions): #Finding the correspondant coordinates of the functional domains in the PDB file structure
            try:
                specific_corresponding_domain_position = Final_Positions[domain_position - 1]
                try:
                    Equivalent_Domain_positions.append(Residues_ID[specific_corresponding_domain_position - 1])
                except:
                    pass
            except:
                pass

        return Equivalent_Domain_positions
def List_of_positions_of_Positive_Sites(file,prob):  # file = Out file from M8 #prob = 95 or 99
        ''' Reads M8 from M7vsM8 Codeml test: Will return the Positions of the selected sites with regards to the ALIGNMENT, later use function to find the equivalent positions in the clean from missing data sequence'''
        length_alignment = []
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
                    if prob == 99 :
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
def Corresponding_positions_missing_notmissing_data(Missing_data,Clean): #Missing_data = WITH missing data ; Clean = NO missing data
    ''' Returns list of the equivalent positions among 2 sequences, in this case the same sequence with and without missing data'''

    alignment_info_global, alignment_global = Global_alignment(Missing_data,Clean)
    Aligned_A, Aligned_B, score, begin, end = alignment_info_global #A = Missing data ('N'); B = Not missing data
    #print(format_alignment(*alignment_global[0]))
    List = equivalent_positions(Missing_data,Clean, Aligned_A, Aligned_B) #All corresponding positions of the sequence with missing data to the sequence without missing data
    return List
def Corresponding_functional_positions(PDB,Full_PDB_sequence,Residues_ID,functional_list):
    alignment_info_global, alignment_global = Global_alignment(PDB,Full_PDB_sequence) #CHECK THIS is in the right order and it should be the other way round
    Aligned_A, Aligned_B, score, begin, end = alignment_info_global
    # print(format_alignment(*alignment_global[0]))
    List = equivalent_positions(PDB,Full_PDB_sequence,Aligned_A,Aligned_B,Residues_ID=Residues_ID,Domain_Positions=functional_list)  # All corresponding functional domain positions in the PDB
    return List
def PROSITE_domains(PDB_sequence):
    "Returns positions in index 1"
    import itertools
    from Bio.ExPASy import ScanProsite
    handle = ScanProsite.scan(seq=PDB_sequence) #can use a sequence (seq=) or a pdbID
    result = ScanProsite.read(handle)
    domains=[]
    for index,segment in enumerate(result,1):
        #if 'score' in segment:
            domains.append(list(range(segment['start'],segment['stop'])))
        #else:
            #pass
    return list(itertools.chain.from_iterable(domains))
def Corresponding_Coordinates_and_labels_PDB_Gene(PDB,Gene,List_of_Positive_Positions,functional_list,Residues_ID,basepath,print_alignment,Gene_name,Clean_positions = None):

    '''Performs the alignments, retrieves the equivalent coordinates and labels among the PDB sequence and the gene for each of the PDB positions.
    Produces a dataframe with different labels for each position, where 'Not'refers to not being positively selected and not functional domain,
    'Selected' stands for positive selected, 'Domain', states that it belongs to the functional domain or 'Selected_and_Domain' which stands for both'''

    #Perform the global alignment
    alignment_info_global, alignment_global = Global_alignment(PDB, Gene)
    Aligned_A, Aligned_B, score, begin, end = alignment_info_global #A = PDB; B = Gene(with or without missing data)
    if print_alignment == 'yes':
        print(format_alignment(*alignment_global[0]))
    else:
        pass

    List_positions = list(range(1, len(Gene) + 1)) #create a list of numbers in range 1 to gene length
    #Extract the corresponding positions of the positive selected sites in the clean of missing data gene sequence
    List_positive_positions = []
    if Clean_positions: #If the gene sequence has been cleaned (removed missing information) and we have the corresponding positions
        for element in List_of_Positive_Positions:
            try:
                List_positive_positions.append(Clean_positions.index(element)) #Append the index of the residue that has the same id as the one in positive selected sites
            except:
                pass
    else: #we don't do anything
        List_positive_positions = List_of_Positive_Positions
    #For the dataframe we can label the positions from 1 to length of gene sequence (with or without missing data)
    positions_dataframe = pd.DataFrame(pd.Series(List_positions))
    positions_dataframe.rename(columns={positions_dataframe.columns[0]: "Gene_Position"}, inplace=True)
    List = equivalent_positions(PDB, Gene, Aligned_A, Aligned_B,Residues_ID = Residues_ID) #List of equivalent positions of each of the gene residues into the PDB structure(residues ID!)
    ###The functional domains require also processing to find their equivalent residues in the PDB

    functional_list = [Residues_ID[index-1] for index in functional_list ] #Index 1 in both, residues ID and Prosite,so fix for python 0 indexing

    #Add the positions to the dataframe
    positions_dataframe['PDB_Position'] = List
    Label_1 = positions_dataframe['Gene_Position'].isin(List_positive_positions) #Check if the position in the gene sequence is + selected
    Label_2 = positions_dataframe['PDB_Position'].isin(functional_list) #Check if the PDB position is a functional domain

    positions_dataframe['Label_1'] = Label_1  # Create a Column were the positive positions have a label True
    positions_dataframe['Label_2'] = Label_2 #Create a column were the functional domain positions recieve the label True
    positions_dataframe['Label_1'] = positions_dataframe['Label_1'].replace(True, 'Selected')
    positions_dataframe['Label_1'] = positions_dataframe['Label_1'].replace(False, 'Not')
    positions_dataframe['Label_2'] = positions_dataframe['Label_2'].replace(True, 'Domain')
    positions_dataframe['Label_2'] = positions_dataframe['Label_2'].replace(False, 'Not')
    positions_dataframe['Label'] = pd.Series(List_positions)
    for index, row in positions_dataframe.iterrows():

        if positions_dataframe.iloc[index, positions_dataframe.columns.get_loc('Label_1')] == 'Selected' and positions_dataframe.iloc[index, positions_dataframe.columns.get_loc('Label_2')] == 'Domain':
            positions_dataframe.iloc[index, positions_dataframe.columns.get_loc('Label')] = 'Selected_and_Domain'
        elif positions_dataframe.iloc[index, positions_dataframe.columns.get_loc('Label_1')] == 'Selected':
            positions_dataframe.iloc[index, positions_dataframe.columns.get_loc('Label')] = 'Selected'
        elif positions_dataframe.iloc[index, positions_dataframe.columns.get_loc('Label_2')] == 'Domain':
            positions_dataframe.iloc[index, positions_dataframe.columns.get_loc('Label')] = 'Domain'
        else:
            positions_dataframe.iloc[index, positions_dataframe.columns.get_loc('Label')] = 'Not'
    positions_dataframe.drop(['Label_1', 'Label_2'], 1, inplace=True)
    Directory = os.path.dirname(basepath)  # Obtain absolute path from the gene file, to create the dataframe there
    Directory=os.path.abspath(os.path.join(os.path.dirname(basepath), '..','Positions_Dataframes'))
    #Gene_name = ntpath.basename(Gene)
    PDB_name = ntpath.basename(basepath)
    positions_dataframe.to_csv(Directory + "/%s_%s_Positions" % (Gene_name,PDB_name), index=False, sep='\t')
    print('Dataframe Ready at %s!' % (Directory))
    return positions_dataframe
def Wrapper_of_all_functions(PDB_sequence,Gene,Chains,M8,List_Domains,Format,prob,missing_data,Residues_ID,basepath,print_alignment,Gene_name):
    '''Calls the function that builds the dataframe of positions that will be read by PyMOL according to according to the optional arguments
    given by the user'''
    List_of_Positive_Positions =List_of_positions_of_Positive_Sites(M8,prob)
    Gene_missing_data= Gene
    #Checking if the user wants to perform the alignment with or without missing data in the gene
    if missing_data == 'no': ###FIX THIS PART!!! The sequence 'Gene' will be translated at this point no matter what
        Clean_protein_sequence = Translate_and_Remove_missing_data(Gene) #Translate and remove missing data from the gene
        Protein_missing_data = Translate_sequence(Gene)  # Translate gene
        Clean_positions = Corresponding_positions_missing_notmissing_data(Protein_missing_data,Clean_protein_sequence) #Find the equivalent positions among the protein with and without missing data
        Dataframe =Corresponding_Coordinates_and_labels_PDB_Gene(PDB_sequence,Clean_protein_sequence, List_of_Positive_Positions,List_Domains,Residues_ID,basepath,print_alignment,Gene_name, Clean_positions)

    else: #Gene_missing_data is our sequence
        Protein_missing_data = Translate_sequence(Gene) #Translate
        Dataframe =Corresponding_Coordinates_and_labels_PDB_Gene(PDB_sequence,Protein_missing_data,List_of_Positive_Positions,List_Domains,Residues_ID,basepath,print_alignment,Gene_name)
        #Corresponding_Coordinates_and_labels_PDB_Gene(PDB,Gene,List_of_Positive_Positions,functional_list,Residues_ID,basepath,print_alignment,Clean_positions = None):
    return Dataframe
def Folders(Genes, folder_name):
    """ Folder for all the generated images It will updated everytime!!! Save the previous folder before running again. Give full path to the Genes sequences file"""
    import os
    import shutil
    basepath = os.path.dirname(Genes)
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
        shutil.rmtree(newpath)  # removes all the subdirectories!
        os.makedirs(newpath,0o777)
def Calling_Pymol():
    """ Call Pymol in each gene and corresponding PDB file in the dataframe"""
    #Read the input fasta file and divide the gene ids and the sequences
    List_protein_names =fasta_to_sequences(Genes,Gene_file_format)[1]
    List_protein_sequences=fasta_to_sequences(Genes,Gene_file_format)[0]
    List_protein_sequences=[Translate_and_Remove_missing_data(sequence) for sequence in List_protein_sequences]
    #Generate the full and filtered summary dataframe with the top suggested homologous proteins from PDB

    PDB_files_dataframe = Download_selected_PDB_files_and_add_resolution_and_percentage_id(Initial_dataframe(List_protein_sequences,List_protein_names))
    #To speed up the process if pdb files already downloaded----> Remember to comment out the create and delete old folder
    #PDB_files_dataframe = pd.read_csv('Full_Blast_results_against_PDB_Filtered.tsv',sep='\t')
    #Read the file with paths to codeml output
    Codeml = pd.read_csv(codeml_output, sep='\s+', header=None)
    Fasta_ids= Codeml.iloc[:,0].tolist()
    M8s= Codeml.iloc[:,1].tolist() #paths
    dictionary=dict(zip(Fasta_ids,M8s))

    for index,row in PDB_files_dataframe.iterrows():

        Gene_name = row["Gene"]
        PDB_file_name =row["PDB_files"].strip("'")
        #Chains =row["Chains"].strip("'") #','.join(list(chains))
        Chains =','.join(list(row['Chains'].strip("'")))#I don't think I want a list here...
        Gene= ''.join([str(fasta.seq) for fasta in SeqIO.parse(open(Genes),Gene_file_format) if fasta.id == Gene_name])#maybe str(fasta.seq)
        if not basepath:
            try:
                PDB_file = os.path.join("PDB_files" + "/pdb%s.ent" % row["PDB_files"].strip("'").lower())
                Residues_ID,PDB_sequence= Extract_sequence_from_PDB(PDB_file,Chains)
            except:
                PDB_file = os.path.join('PDB_files'+ "/%s.cif" % row["PDB_files"].strip("'").lower())
                Residues_ID, PDB_sequence = Extract_sequence_from_PDB(PDB_file, Chains)
        else:
            try:
                PDB_file = os.path.join(basepath + "/PDB_files" + "/pdb%s.ent" % row["PDB_files"].strip("'").lower())
                Residues_ID,PDB_sequence= Extract_sequence_from_PDB(PDB_file,Chains)
            except:
                PDB_file = os.path.join(basepath + '/PDB_files'+ "/%s.cif" % row["PDB_files"].strip("'").lower())
                Residues_ID, PDB_sequence = Extract_sequence_from_PDB(PDB_file, Chains)
        #Extract positions that belong to domain in index 1

        M8 = ''.join([path for id,path in dictionary.items() if id==Gene_name])
        try:
            List_domains = PROSITE_domains(PDB_sequence)
            Data = Wrapper_of_all_functions(PDB_sequence, Gene, Chains, M8, List_domains, Gene_file_format, prob,missing_data, Residues_ID, PDB_file, print_alignment,Gene_name)
        except:
            pass
            #List_domains = []
            #Data = Wrapper_of_all_functions(PDB_sequence, Gene, Chains, M8, List_domains, Gene_file_format, prob,missing_data, Residues_ID, PDB_file, print_alignment,Gene_name)

if __name__ == "__main__":

    #Folders(Genes,"PDB_files")
    Folders(Genes,'Positions_Dataframes')
    Folders(Genes,'LYS_Pymol_Images')
    Calling_Pymol()
