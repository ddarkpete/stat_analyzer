from Bio.PDB import *
from Bio.PDB.Atom import Atom 
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
from Bio.PDB.PDBIO import PDBIO
from collections import OrderedDict

import os
from os.path import join

def percentage(x , _all):
    return (float(x) / float(_all)) * 100


def stat_reader():
    file = "stats.txt"
    stats = {}
    with open(file , 'r') as input_stats:
        for line in input_stats:
            splited =  line.split('\t')
            splited[2] = splited[2].rstrip()
            temp_tulp = (splited[0] , splited[1])
            stats[temp_tulp] = float(splited[2])
    return stats


path = "./models_to_score/"

#model_stats = {}

standard_aa_names = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS",
                     "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL",
                     "TRP", "TYR"]

stats = stat_reader()
pair_counter = 0

all_models_stats = []

files = [f for f in os.listdir(path) if os.path.isfile(join(path,f)) and '.pdb' in f]

for pdb_file in files:#os.listdir(path):
    
    print pdb_file + '\n'
    parser = PDBParser()
    struct = parser.get_structure('structure', path + pdb_file)
    chains = list(struct.get_chains())
    compares = []#storing compares that are already done
    #analyzed_count += 1
    model_stats = {}
    
    for ch1 in range(0,len(chains)):
        for ch2 in range(ch1 + 1,len(chains)):
            checklist = [ chains[ch1].get_full_id()[2] , chains[ch2].get_full_id()[2] ]
            checklist2 =  [ chains[ch2].get_full_id()[2], chains[ch1].get_full_id()[2] ]
            if chains[ch1].get_full_id()[2] != chains[ch2].get_full_id()[2] and not checklist in compares and not checklist2 in compares:

                comparsion = [ chains[ch1].get_full_id()[2] , chains[ch2].get_full_id()[2] ]
                compares.append(comparsion)#appending comprasion to already done comprasions

                chain1_atms = list(chains[ch1].get_atoms())
                chain2_atms = list(chains[ch2].get_atoms())
                
                for atm1 in chain1_atms:
                    if atm1.get_name() == 'CB' and atm1.get_parent().get_resname() in standard_aa_names:
                        #print 'CB atom 1'
                        for atm2 in chain2_atms:
                            if atm2.get_name() == 'CB' and atm2.get_parent().get_resname() in standard_aa_names:
                                #print 'cb atom 2'
                                if atm1 - atm2 <= 10.0:
                                    pair1 = (atm1.get_parent().get_resname() , atm2.get_parent().get_resname()) 
                                    pair2 = (atm2.get_parent().get_resname() , atm1.get_parent().get_resname())

                                    if pair1 in stats:
                                        if pair1 in model_stats:
                                            model_stats[pair1] += 1
                                            pair_counter += 1
                                        else:
                                            model_stats[pair1] = 1
                                            pair_counter += 1
                                    elif pair2 in stats:
                                        if pair2 in model_stats:
                                            model_stats[pair2] += 1
                                            pair_counter += 1
                                        else:
                                            model_stats[pair2] = 1
                                            pair_counter += 1
                                    else:
                                        model_stats[pair1] = 1
                                        pair_counter += 1
    for pair in model_stats:
    	model_stats[pair] = percentage(model_stats[pair] , pair_counter)
    	if model_stats[pair] in stats:
    		model_stats[pair] = model_stats[pair] - stats[pair]
    for pair in stats:
    	if pair not in model_stats:
    		model_stats[pair] = 0 - stats[pair]

    	

    all_models_stats.append(model_stats)

#TODO
# - Struktura dla par - licznosc modeli 
# - Odrazu porownanie i print roznic licznosci par / zapis do pliku ?
# - 

for am in all_models_stats:
    for ms , ms_val in am.iteritems():
    	print str(ms) +  ' ' + str(ms_val)
    print ""




