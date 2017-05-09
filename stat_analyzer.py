from Bio.PDB import *
from Bio.PDB.Atom import Atom 
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
from Bio.PDB.PDBIO import PDBIO
from collections import OrderedDict

#import chimera
import os

path = ""

#parser = PDBParser()
files = [f for f in os.listdir('.') if os.path.isfile(f) and '.pdb' in f]
resis_on_close = []
res_ids_onclose = {}
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
                chain1_resis = []#lists of residues involved in contact between two chains in chain_1
                chain2_resis = []#same as above
                chain1_close_atms = []
                chain2_close_atms = []

                for atm1 in chain1_atms:
                    for atm2 in chain2_atms:
                        if atm1 - atm2 <= 10.0:
                            #print atm1.get_name() + ' ' + atm2.get_name()
                            if atm1.get_parent().get_resname() not in chain1_resis:
                                chain1_resis.append(atm1.get_parent().get_resname())
                            
                            if atm2.get_parent().get_resname() not in chain2_resis:
                                chain2_resis.append(atm2.get_parent().get_resname())
                            
                            if atm1.get_name() not in chain1_close_atms:
                                chain1_close_atms.append(atm1.get_name())
                            
                            if atm2.get_name() not in chain2_close_atms:
                                chain2_close_atms.append(atm2.get_name()) 
                
                if len(chain1_resis) == 0 or len(chain2_resis) == 0:
                    print ' No atoms in distance <= 10 Angstremes for ' + chains[ch1].get_full_id()[2] + ' and ' + chains[ch2].get_full_id()[2] + ' chains' + '\n'
                else:
                    print 'Residues and atoms in distance <= 10 Angstremes for ' + chains[ch1].get_full_id()[2] + ' and ' + chains[ch2].get_full_id()[2] + ' chains' + '\n'
                    print chain1_resis #close residues from 1st chain
                    print ""
                    print chain1_close_atms #close atoms from 1st chain
                    print "\n"
                    print chain2_resis
                    print ''
                    print chain2_close_atms

                    
                    for r1 in chain1_resis:
                        if r1 not in res_ids_onclose:
                            res_ids_onclose[r1] = 1
                        else:
                            res_ids_onclose[r1] += 1
                    #print ''
                    #print chain2_close
                    for r2 in chain2_resis:
                        if r2 not in res_ids_onclose:
                            res_ids_onclose[r2] = 1
                        else:
                            res_ids_onclose[r2] += 1
                    #print ''
                    break
            

sorted_by_count = OrderedDict(sorted(res_ids_onclose.items(), key=lambda x: x[1]))
#sorted_by_count = reverse(sorted_by_count)

print sorted_by_count # slownik reszt wraz z liczba ich wystepowan 