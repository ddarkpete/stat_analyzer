'''
import warnings
warnings.filterwarnings("ignore")
'''
from Bio.PDB import *
from Bio.PDB.Atom import Atom 
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
from Bio.PDB.PDBIO import PDBIO
from collections import OrderedDict
from math import log10

import os
from os.path import join

class Model:
    def __init__(self):
        self.name = ""
        self.stats = {}
        self.score = 0

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

def aa_observed_reader():
    file = "aa_observed.txt"
    with open(file , 'r') as aas:
        buff = aas.read()
        observed_aas = eval(buff)
    return observed_aas


path = "./models_to_score/test/test/"

#model.stats = {}

standard_aa_names = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS",
                     "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL",
                     "TRP", "TYR"]


print str(len(standard_aa_names)) + 'aminos'

stats = stat_reader()
print stats
aa_observed = aa_observed_reader()
print aa_observed
print'************************************'
pair_counter = 0.0

all_models_stats = []

files = [f for f in os.listdir(path) if os.path.isfile(join(path,f)) and '.pdb' in f]

for pdb_file in files:#os.listdir(path):
    
    print pdb_file + '\n'
    parser = PDBParser()
    struct = parser.get_structure('structure', path + pdb_file)
    chains = list(struct.get_chains())
    #print str(len(chains))
    compares = []#storing compares that are already done
    #analyzed_count += 1
    model = Model()
    model.name = pdb_file
    
    for ch1 in range(0,len(chains)):
        for ch2 in range(ch1 + 1,len(chains)):
            checklist = [ chains[ch1].get_full_id()[2] , chains[ch2].get_full_id()[2] ]
            checklist2 =  [ chains[ch2].get_full_id()[2], chains[ch1].get_full_id()[2] ]
            if chains[ch1].get_full_id()[2] != chains[ch2].get_full_id()[2] and not checklist in compares and not checklist2 in compares:

                comparsion = [ chains[ch1].get_full_id()[2] , chains[ch2].get_full_id()[2] ]
                compares.append(comparsion)#appending comprasion to already done comprasions

                chain1_atms = list(chains[ch1].get_atoms())
                chain2_atms = list(chains[ch2].get_atoms())
                #print str(len(chain1_atms)) + ' ' + str(len(chain2_atms))
                
                for atm1 in chain1_atms:
                    if atm1.get_name() == 'CB' and atm1.get_parent().get_resname() in standard_aa_names:
                        #print 'CB atom 1'
                        for atm2 in chain2_atms:
                            if atm2.get_name() == 'CB' and atm2.get_parent().get_resname() in standard_aa_names:
                                #print 'cb atom 2'
                                if atm1 - atm2 <= 20.0:
                                    #print 'in 10 angs'
                                    pair1 = (atm1.get_parent().get_resname() , atm2.get_parent().get_resname()) 
                                    pair2 = (atm2.get_parent().get_resname() , atm1.get_parent().get_resname())

                                    if pair1 in stats:
                                        if pair1 in model.stats:
                                            model.stats[pair1] += 1
                                            pair_counter += 1
                                        else:
                                            model.stats[pair1] = 1
                                            pair_counter += 1
                                    elif pair2 in stats:
                                        if pair2 in model.stats:
                                            model.stats[pair2] += 1
                                            pair_counter += 1
                                        else:
                                            model.stats[pair2] = 1
                                            pair_counter += 1
    blacklist = []
    for pair in model.stats:
        #model.stats[pair] = percentage(model.stats[pair] , pair_counter)
        if pair in stats:#a co jak nie ?
            model.stats[pair] = log10(float(model.stats[pair]) / float(stats[pair]))
        else:
            #print pair
            blacklist.append(pair)
    
    for b in blacklist:
        del model.stats[b]

    model.score = sum(model.stats.values())
    '''
    for pair in stats:
        if pair not in model.stats:
            model.stats[pair] = 0 - stats[pair]#no i co w tym wypadku ?????
    '''
        

    all_models_stats.append(model)

#TODO
# - Struktura dla par - licznosc modeli 
# - Odrazu porownanie i print roznic licznosci par / zapis do pliku ?
# - 

for am in all_models_stats:
    print "\n" + am.name + "\n"
    for ms , ms_val in am.stats.iteritems():
        print str(ms) +  ' ' + str(ms_val)
    print "\n model score : " + str(am.score) + '\n' 

print str(len(all_models_stats))




