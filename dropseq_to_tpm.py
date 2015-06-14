"""
File to convert the dropseq count files into TPM.
"""
from __future__ import division, print_function
import sys;
import numpy as np;
import re;

if(len(sys.argv) == 1 or sys.argv[1] == '-h'):
    #print help and exit
    print("Usage:  python dropseq_to_tpm.py <count matrix filename>")
    sys.exit();

filename = sys.argv[1];
gene_labels = list();
with open(filename, 'r') as fin:
    firstline = fin.readline();
    for line in fin:
        gene_label = line.split('\t')[0];
        gene_labels.append(gene_label);

col_labels = firstline.strip().split('\t');
num_cols = len(col_labels);

count_matrix = np.loadtxt(filename, dtype=np.float32,skiprows=1,usecols=np.arange(1,num_cols));
col_labels = col_labels[1:];  #First column is just labeled "gene"

start_regex = re.compile(r':\d+\-');
end_regex = re.compile(r'\-\d+:');

transcript_starts = np.zeros(len(gene_labels));
transcript_ends = np.zeros(len(gene_labels));
gene_names = list();

for i, gene_label in enumerate(gene_labels):
    start = start_regex.search(gene_label).group()[1:-1];
    transcript_starts[i] = start;
    end = end_regex.search(gene_label).group()[1:-1];
    transcript_ends[i] = end;
    gene_name = gene_label[gene_label.rfind(':')+1:];
    gene_names.append(gene_name);

transcript_lengths = transcript_ends - transcript_starts;
transcript_lengths.shape = (transcript_lengths.size, 1);

#Use transcript_lengths and counts to compute tpm
total_reads = np.sum(count_matrix, axis=0, keepdims=True);

#convert count_matrix to TPM in place
tpm = count_matrix;
np.divide(tpm, total_reads, tpm);
np.divide(tpm, transcript_lengths, tpm);

norm_factor = np.sum(tpm, axis=0, keepdims=True);
np.divide(tpm, norm_factor, tpm);
np.multiply(tpm, 1e6, tpm);

#convert to log(tpm+1)
np.add(tpm, 1, tpm);
np.log(tpm, tpm);

#Output the results

#Determine output file name
i = filename.rfind('.');
if(i != -1):
    outfile = filename[0:i]+".tpm";
else:
    outfile = filename + ".tpm";

with open(outfile,'w') as fout:
    fout.write("\t".join(col_labels) + "\n");
    for i,gene in enumerate(gene_names):
        tpm_str = ["{:.4g}".format(x) for x in tpm[i,:]];
        fout.write(gene + '\t' + '\t'.join(tpm_str) + '\n');
        print(i);


