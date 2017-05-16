from Bio.PDB import *
from Bio.PDB.Atom import Atom 
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
from Bio.PDB.PDBIO import PDBIO
from collections import OrderedDict

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

stats = stat_reader()
for pair , stat in stats.iteritems():
	print str(pair) + '    ' + str(stat)
