# scikit-learn bootstrap

from collections import Counter
from itertools import combinations, chain
from warnings import warn
from math import isinf, inf
import itertools
from itertools import permutations
import matplotlib.pyplot as plt
import re 
import itertools
import sys
import ast



def genus_comb(genus_list):
	"""
    Generates all possible combinations of the species list provided to script.
    Args:
        genus_list: list of microbial xml files given as an input 
    Returns:
        float: a list containing all combinations of species (no repetition) 
    """

	powerSet = []
	n = len(genus_list)+1
	for k in range(n):
		powerSet.extend(itertools.combinations(genus_list, k))

	return powerSet


def main(argv):

	genus_string = argv
	genus_trial = genus_string.split(",")  
	#print(genus_trial)
	genus_comb_list = genus_comb(genus_trial)
	genus_comb_list = [', '.join(map(str, x)) for x in genus_comb_list]
	command_len = len(genus_comb_list) - len(genus_trial)
	mycommandlist = [None]*(command_len)
	num = 0 
	for i in range(len(genus_trial)+1,len(genus_comb_list)):
		
		comb = genus_comb_list[i] 
		
		mycommandlist[num] = comb
		print(mycommandlist[num])
		num = num +1 
	


if __name__ == '__main__':
	if len(sys.argv) > 1:
		main(sys.argv[1])



