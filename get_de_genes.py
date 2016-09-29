#!/usr/bin/env python
'''
Created on 1 Apr 2016
@author: steve
'''
"""
A class for reference sequences  - stores header:seq pairs in a dictionary
{header:seq}
"""
import time
import argparse
import sys

class Ref_Seq(object):
    def __init__(self):
        self._internal_dict = {}
    
    def __setitem__(self, header, sequence):
        self._internal_dict[header]=sequence
    
    def __getitem__(self, header):
        return self._internal_dict[header]
    
    def __iter__(self):
        return self._internal_dict.iteritems()
    
    def __len__(self):
        return len(self._internal_dict)

    def __contains__(self, header):
    	return header in self._internal_dict


    def headers(self):
        return self._internal_dict.keys()

    def sequences(self):
        return self._internal_dict.values()    
    

    def load_ref_file(self, ref_file):
        """
        Load ref file in fasta format
    
        Product ref_dict --> header:sequence
        """

        start = time.clock()
        #ref_dict = {}
        ref_count = 0
        loaded_ref = open(ref_file, 'rU')
        full_len_seq = ''
        key = ''
        first_header = True
        for line in loaded_ref:
            if line[0] == '>' and full_len_seq == '':
                key = line.strip()
                ref_count += 1
                if first_header:
                    first_header = False
            elif line[0] == '>' and full_len_seq != '':
                self._internal_dict.update({key: full_len_seq})
                key = line.strip()
                full_len_seq = ''
                ref_count += 1
            elif line[0] == '' and full_len_seq != '':
                self._internal_dict.update({key: full_len_seq})
                key = line.strip()
                full_len_seq = ''
            elif line[0] == '':
                pass
            else:
                full_len_seq += line.strip().upper()
    
        self._internal_dict.update({key: full_len_seq})
    
        print '\n----{0} reference sequences loaded for alignment----'\
            .format(ref_count)
        if len(self._internal_dict) ==1:
            print "\n{0} length = {1} bp".format(ref_file.split('/')[-1],
                                             len(full_len_seq))
        print "\nReference sequence loading time = "\
+ str((time.clock() - start)) + " seconds\n"

def parse_counts_matrix(in_file):
	"""
	Get gene names from counts matrix
	"""
	with open(in_file,'r') as f:
		gene_list=[]
		for line in f:
			gene_list.append(">"+line.strip().split('\t')[0])
	return gene_list

def main():
	parser = argparse.ArgumentParser( description='Write de genes to FASTA file')
	## output file to be written
	parser.add_argument('-a', '--assembly', type=str, required=True, help='Path to unigene assembly file')
	parser.add_argument('-d', '--de', type=str, required=True, help='Path to de file')
	parser.add_argument('-o', '--output', type=str, required=False, help='Output file to be created.  Default = STDOUT')
	args = parser.parse_args()
	fout = sys.stdout
	if args.output is not None:
		fout = open(args.output, 'wt')
	gene_list = parse_counts_matrix(args.de)
	assembly=Ref_Seq()
	assembly.load_ref_file(args.assembly)
	de_gene_dict={}
	for gene in gene_list:
		if gene in assembly:
			seq=assembly[gene]
			de_gene_dict[gene]=seq
	for key, value in de_gene_dict.iteritems():
		fout.write(key+"\n")
		lines = (value[n:n + 64] for n in range(0, len(value), 64))
		for i in lines:
			fout.write(i+"\n")
	fout.close()


if __name__ == '__main__':
	main()

