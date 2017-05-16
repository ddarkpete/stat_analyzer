import os
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


path = ""

#parser = PDBParser()
files = [f for f in os.listdir('.') if os.path.isfile(f) and '.pdb' in f]
resis_on_close = []
comprasions = []
res_ids_onclose = {}
pairs_occurence = {}
analyzed_count = 0



for pdb_file in files:#os.listdir(path):
    
    print pdb_file + '\n'
    parser = PDBParser()
    struct = parser.get_structure('structure', pdb_file)
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
                
                '''
                chain1_resis = []#lists of residues involved in contact between two chains in chain_1
                chain2_resis = []#same as above
                chain1_close_atms = []
                chain2_close_atms = []
                '''

                for atm1 in chain1_atms:
                    if atm1.get_name() == 'CB':
                        #print 'CB atom 1'
                        for atm2 in chain2_atms:
                            if atm2.get_name() == 'CB':
                                #print 'cb atom 2'
                                if atm1 - atm2 <= 10.0:
                                    #print atm1.get_name() + ' ' + atm2.get_name()
                                    if atm1.get_parent().get_resname() not in comp.chain1_resis:
                                        comp.chain1_resis.append(atm1.get_parent().get_resname())
                                    
                                    if atm2.get_parent().get_resname() not in comp.chain2_resis:
                                        comp.chain2_resis.append(atm2.get_parent().get_resname())
                                    
                                    if atm1.get_name() not in comp.chain1_close_atms:
                                        comp.chain1_close_atms.append(atm1.get_name())
                                    
                                    if atm2.get_name() not in comp.chain2_close_atms:
                                        comp.chain2_close_atms.append(atm2.get_name())
                                    
                                    pair1 = (atm1.get_parent().get_resname() , atm2.get_parent().get_resname()) 
                                    pair2 = (atm2.get_parent().get_resname() , atm1.get_parent().get_resname())

                                    if pair1 in pairs_occurence:
                                        pairs_occurence[pair1] += 1
                                    elif pair2 in pairs_occurence:
                                        pairs_occurence[pair2] += 1
                                    else:
                                        pairs_occurence[pair1] = 1
                
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
#print pairs_occurence
for data in pairs_occurence:
    print str(data[0]) + ' occures ' + str(data[1]) + ' times'