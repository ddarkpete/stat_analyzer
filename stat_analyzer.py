import os
from os.path import join
#import chimera

from Bio.PDB import *
from Bio.PDB.Atom import Atom 
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
from Bio.PDB.PDBIO import PDBIO
from collections import OrderedDict



class Comprassion:
    def __init__(self):
        self.chain1_resis = []
        self.chain2_resis = []
        self.chain1_close_atms = []
        self.chain2_close_atms = []

def percentage(x , _all):
    return (float(x) / float(_all)) * 100

def convert_to_percentage(occ_list, _all):
    save_file = open('stats.txt', 'w')
    for data in occ_list:
        list(data)
        #print data
        data[1] = percentage(data[1], _all)
        save_file.write(data[0][0] + '\t' + data[0][1] +  '\t' + str(data[1]) + '\n')
    

standard_aa_names = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS',
                     'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL',
                     'TRP', 'TYR']

aa_observed = {'ALA' : 0 , 'CYS' : 0 , 'ASP' : 0 , 'GLU' : 0 , 'PHE' : 0 , 'GLY' : 0 , 'HIS' : 0 , 'ILE' : 0 , 'LYS' : 0 ,
                     'LEU' : 0 , 'MET' : 0 , 'ASN' : 0 , 'PRO' : 0 , 'GLN' : 0 , 'ARG' : 0 , 'SER' : 0 , 'THR' : 0 , 'VAL' : 0 ,
                     'TRP' : 0 , 'TYR' : 0 }


path = "./pdbs/"

#parser = PDBParser()
print os.listdir(path)
files = [f for f in os.listdir(path) if os.path.isfile(join(path,f)) and '.pdb' in f]
print files
resis_on_close = []
comprasions = []
res_ids_onclose = {}
pairs_occurence = {}
analyzed_count = 0
pair_counter = 0



for pdb_file in files:#os.listdir(path):
    
    print pdb_file + '\n'
    parser = PDBParser()
    struct = parser.get_structure('structure', path + pdb_file)
    chains = list(struct.get_chains())
    compares = []#storing compares that are already done
    analyzed_count += 1
    
    for ch1 in range(0,len(chains)):
        for ch2 in range(ch1 + 1,len(chains)):
            checklist = [ chains[ch1].get_full_id()[2] , chains[ch2].get_full_id()[2] ]
            checklist2 =  [ chains[ch2].get_full_id()[2], chains[ch1].get_full_id()[2] ]
            if chains[ch1].get_full_id()[2] != chains[ch2].get_full_id()[2] and not checklist in compares and not checklist2 in compares:

                comparsion = [ chains[ch1].get_full_id()[2] , chains[ch2].get_full_id()[2] ]
                compares.append(comparsion)#appending comprasion to already done comprasions

                chain1_atms = list(chains[ch1].get_atoms())
                chain2_atms = list(chains[ch2].get_atoms())

                comp = Comprassion()
                
                

                for atm1 in chain1_atms:
                    if atm1.get_name() == 'CB' and atm1.get_parent().get_resname() in standard_aa_names:
                        #print 'CB atom 1'
                        for atm2 in chain2_atms:
                            if atm2.get_name() == 'CB' and atm2.get_parent().get_resname() in standard_aa_names:
                                #print 'cb atom 2'
                                if atm1 - atm2 <= 10.0:
                                    #print atm1.get_name() + ' ' + atm2.get_name()
                                    if atm1.get_parent().get_resname() not in comp.chain1_resis:
                                        comp.chain1_resis.append(atm1.get_parent().get_resname())
                                        aa_observed[atm1.get_parent().get_resname()]+=1
                                    if atm2.get_parent().get_resname() not in comp.chain2_resis:
                                        comp.chain2_resis.append(atm2.get_parent().get_resname())
                                        aa_observed[atm2.get_parent().get_resname()]+=1
                                    if atm1.get_name() not in comp.chain1_close_atms:
                                        comp.chain1_close_atms.append(atm1.get_name())
                                    
                                    if atm2.get_name() not in comp.chain2_close_atms:
                                        comp.chain2_close_atms.append(atm2.get_name())
                                    
                                    pair1 = (atm1.get_parent().get_resname() , atm2.get_parent().get_resname()) 
                                    pair2 = (atm2.get_parent().get_resname() , atm1.get_parent().get_resname())

                                    if pair1 in pairs_occurence:
                                        pairs_occurence[pair1] += 1
                                        pair_counter += 1
                                    elif pair2 in pairs_occurence:
                                        pairs_occurence[pair2] += 1
                                        pair_counter += 1
                                    else:
                                        pairs_occurence[pair1] = 1
                                        pair_counter += 1
                
                if len(comp.chain1_resis) == 0 or len(comp.chain2_resis) == 0:
                    print ' No atoms in distance <= 10 Angstremes for ' + chains[ch1].get_full_id()[2] + ' and ' + chains[ch2].get_full_id()[2] + ' chains' + '\n'
                else:
                    
                    comprasions.append(comp)

                    print 'Residues and atoms in distance <= 10 Angstremes for ' + chains[ch1].get_full_id()[2] + ' and ' + chains[ch2].get_full_id()[2] + ' chains' + '\n'
                    print comp.chain1_resis #close residues from 1st chain
                    print ""
                    print comp.chain1_close_atms #close atoms from 1st chain
                    print "\n"
                    print comp.chain2_resis
                    print ''
                    print comp.chain2_close_atms

                    
                    for r1 in comp.chain1_resis:
                        if r1 not in res_ids_onclose:
                            res_ids_onclose[r1] = 1
                        else:
                            res_ids_onclose[r1] += 1
                    #print ''
                    #print chain2_close
                    for r2 in comp.chain2_resis:
                        if r2 not in res_ids_onclose:
                            res_ids_onclose[r2] = 1
                        else:
                            res_ids_onclose[r2] += 1
                    #print ''
                    break
  
print'\n\n\n\n'          

#sorted_by_count = OrderedDict(sorted(res_ids_onclose.items(), key=lambda x: x[1]))
#sorted_by_count = reverse(sorted_by_count)
pairs_occurence = sorted(pairs_occurence.items(), key=lambda kv: kv[1], reverse=True)

#print sorted_by_count # slownik reszt wraz z liczba ich wystepowan 
print'\n\n\n\n'
print str(pair_counter)
print'\n\n\n\n'
sumka = 0
pairs_occurence = [list(elem) for elem in pairs_occurence]
for data in pairs_occurence:
    print str(data[0][0]) + ' - ' + str(data[0][1]) + ' occures ' + str(data[1]) + ' times.' + 'It is ' + str(percentage(data[1],pair_counter)) + ' % of all pairs'
    sumka += data[1]
#print pairs_occurence_inlist
print str(sumka)
convert_to_percentage(pairs_occurence,pair_counter)

''' TODO

-Do przekazania słownik do drugiego pliku z licznośćią aminokwasów
-dalsza implementacja wzoru w sensie tego wyliczjacego expected
-poten ten logarytmem