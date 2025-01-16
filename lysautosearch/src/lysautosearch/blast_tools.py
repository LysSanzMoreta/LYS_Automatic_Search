import pandas as pd
from urllib import request
import re
import warnings
import ast
import os


def blast_dataframe(protein_sequences_missing_data:list,protein_sequences_no_missing_data:list,protein_id_names:list,results_dir:str):
    """Builds a dataframe with all the necessary information from the blast results"""
    dictionary ={}
    for protein_missing_data,protein_nomissing_data,id in zip(protein_sequences_missing_data,protein_sequences_no_missing_data,protein_id_names):
        PDBs_List,chains,percentID,sequences_missing_data,sequences_no_missing_data = blast_extract_PDB(protein_missing_data.replace("*",""),protein_nomissing_data.replace("*",""))
        if not PDBs_List:#to avoid empty information
            pass
        else:
            dictionary[id] = PDBs_List,chains,percentID,sequences_missing_data,sequences_no_missing_data
    df = pd.DataFrame.from_dict(dictionary,orient='index')
    df.index.name = 'Genes'
    df.reset_index(inplace=True)
    df.columns = ['Gene','PDB_files','Chains','PercentID',"Prot_seq_clean","Prot_seq_full"]
    df.to_csv("{}/Full_Blast_results_against_PDB.tsv".format(results_dir),sep='\t',index=False)
    return df


def blast_extract_PDB(protein_no_missing_data:str,protein_missing_data:str): #Translated gene sequence
    """Using the protein sequences, blast them against the PDB database and parse the result, only significant e-values: 1e-6

    perhaps useful: https://github.com/rcsb/py-rcsb-api
    reference url:   #https://search.rcsb.org/query-editor.html?json=%7B%22query%22%3A%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22sequence%22%2C%22parameters%22%3A%7B%22evalue_cutoff%22%3A1%2C%22identity_cutoff%22%3A0.9%2C%22sequence_type%22%3A%22protein%22%2C%22value%22%3A%22MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLPARTVETRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPGCMNCKCVIS%22%7D%7D%2C%22request_options%22%3A%7B%22scoring_strategy%22%3A%22sequence%22%7D%2C%22return_type%22%3A%22polymer_instance%22%7D


    JSON query format:
        {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": 1,
                "identity_cutoff": 0.9,
                "sequence_type": "protein",
                #"value": "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLPARTVETRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPGCMNCKCVIS",
                "value": f"{sequence}"
            }
        },
        "request_options": {
            "scoring_strategy": "sequence"
        },
        "return_type": "polymer_instance"
        }
        ######################################

        "return_type": ["structures","polymer_entities","assemblies","non_polymer_entities","molecular_definitions","polymer_instance]
    """


    PDB_ids=[]
    chains=[]
    identities=[]
    sequences_no_missing_data = []
    sequences_missing_data = []


    evalue_cutoff = 0.1 #lower values (closer to zero, more significant match)
    identity_cutoff = 0 #need to keep low for some reason
    #working_seq = "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLPARTVETRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPGCMNCKCVIS"
    #sequence = "MSNLRDALLNGLKNLPGTREFHLHVLVSSPRKHGGLYPFSXXTSMLYLQDVLILLSEQSTPDSPRIFVAAIEACVYNVAATSSAILYISKVDSTGQASYPSPTATLVKSFLFYYADPATRPIDVDHLWIQLFARAQKQYVFPNSSEYEGKKWLNDMKLCSWWKRIYTDVSTELASRLPSTALVKLYYLIPGYNTTXAEQSLKIASSSTGNDNQNAPMWIYGHPYSQTDIPLPCPQDTGSSSSIVRNLGHYIPSFDDDPKSRFMDEIAYTTEKEGIKSPPRKRARTISSSRHPDAVPVSGSAAAEDDEAEKPGSQDRPLGELSKVSPDEFWERMSFRQECVSGAVTGFFTMGISFKAADEDGSRSTVSPLAPQAGQVSSNVNKRVMTSLMTGVEFSTKERAINATESIENAIKGLCEGIALVPVPAIQPTRPKSTSNQEXXXRQSSDPEVXSTYLAPPRTPPRRNRELLPDISPNPFPEPEPTLETYNSHIYGSVLAVRKKKKRAD"
    url = f"https://search.rcsb.org/rcsbsearch/v2/query?json=%7B%22query%22%3A%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22sequence%22%2C%22parameters%22%3A%7B%22evalue_cutoff%22%3A{evalue_cutoff}%2C%22identity_cutoff%22%3A{identity_cutoff}%2C%22sequence_type%22%3A%22protein%22%2C%22value%22%3A%22{protein_missing_data}%22%7D%7D%2C%22request_options%22%3A%7B%22scoring_strategy%22%3A%22sequence%22%7D%2C%22return_type%22%3A%22polymer_instance%22%7D"
    data = request.urlopen(url).read()
    if data:
        data = data.decode("utf-8").replace("b","").replace("\n","")
        data = ast.literal_eval(data) #convert to dict
        for result in data["result_set"]:
            pdb_id,chain = result["identifier"].split(".")
            score = result["score"] # i am not sure what is this score
            PDB_ids.append(pdb_id)
            chains.append(chain)
            identities.append(score)
            sequences_no_missing_data.append(str(protein_no_missing_data))
            sequences_missing_data.append(str(protein_missing_data))

    else:
        warnings.warn("Sequence not found, skipping")

    return PDB_ids,chains, identities,sequences_no_missing_data,sequences_missing_data