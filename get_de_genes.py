#!/usr/bin/env python

"""
Using a edgeR matrix file of DE genes, extracts these genes from unigene transcriptome, and can generate length and
factor labelling files if required.
"""

import time
import argparse
import sys


class RefSeq(object):
    def __init__(self):
        self._internal_dict = {}

    def __setitem__(self, header, sequence):
        self._internal_dict[header] = sequence

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
        # ref_dict = {}
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

        print '\n{0} transcripts in transciptome\n' \
            .format(ref_count)



def parse_counts_matrix(in_file):
    """
	Get gene names from counts matrix
	"""
    with open(in_file, 'r') as f:
        gene_list = []
        line_count = 0
        for line in f:
            if line_count > 0:
                gene_list.append(">" + line.strip().split('\t')[0])
            line_count += 1
    print "\n{0} differentially expressed transcripts in matrix file".format(line_count - 1)
    return gene_list


def write_factor_labeling_file(gene_list, out_put):
    with open(out_put, 'w') as f:
        for gene in gene_list:
            f.write("{0}\t{1}\n".format(out_put.split('.')[0], gene[1:]))
    f.close()


def write_gene_lengths(args, assembly):
    fout = open(args.output + '_lengths.txt', 'wt')
    fout.write('gene\tlength\n')
    len_count = 0
    gene_count = 0
    for key, value in assembly:
        len_count += (len(value))
        gene_count += 1
        fout.write("{0}\t{1}\n".format(key[1:], len(value)))
    print "Average transcript length = {0}bp\n".format(len_count / gene_count)


def write_de_gene_fasta(assembly, fout, gene_list):
    de_gene_dict = {}
    for gene in gene_list:
        if gene in assembly:
            seq = assembly[gene]
            de_gene_dict[gene] = seq
    for key, value in de_gene_dict.iteritems():
        fout.write(key + "\n")
        lines = (value[n:n + 64] for n in range(0, len(value), 64))
        for i in lines:
            fout.write(i + "\n")
    fout.close()


def main():

    print " _____      _    ______ _____   _____                      "
    print "|  __ \    | |   |  _  \  ___| |  __ \                     "
    print "| |  \/ ___| |_  | | | | |__   | |  \/ ___ _ __   ___  ___ "
    print "| | __ / _ \ __| | | | |  __|  | | __ / _ \ '_ \ / _ \/ __|"
    print "| |_\ \  __/ |_  | |/ /| |___  | |_\ \  __/ | | |  __/\__ \ "
    print " \____/\___|\__| |___/ \____/   \____/\___|_| |_|\___||___/\n"
    parser = argparse.ArgumentParser(description='Write de genes to FASTA file')
    ## output file to be written
    parser.add_argument('-a', '--assembly', type=str, required=True, help='Path to unigene assembly file')
    parser.add_argument('-d', '--de', type=str, required=True, help='Path to de file')
    parser.add_argument('-go', '--go', action='store_true',
                        help='Generate gene lengths and factor labelling files for GO enrichment analysis')
    parser.add_argument('-o', '--output', type=str, required=False, help='Output file to be created.  Default = STDOUT')
    args = parser.parse_args()
    fout = sys.stdout
    gene_list = parse_counts_matrix(args.de)
    if args.output is not None:
        fout = open(args.output + '.fa', 'wt')
        if args.go:
            write_factor_labeling_file(gene_list, args.output + '_factor.txt')

    assembly = RefSeq()
    assembly.load_ref_file(args.assembly)
    write_de_gene_fasta(assembly, fout, gene_list)
    if args.go:
        write_gene_lengths(args, assembly)


if __name__ == '__main__':
    main()
