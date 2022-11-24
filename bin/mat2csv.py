import csv
import gzip
import os
import scipy.io
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Parser for Mat2Csv')
parser.add_argument('--hto_mat_folder', help= 'cellranger output folder which contains HTO (antibody tag) count matrix in mtx format and other meta data.')
parser.add_argument('--gz', help= 'Whether files are in .gz format or not', default=True)

args = parser.parse_args()
hto_mat_folder = args.hto_mat_folder

hto_mat_mtx = scipy.io.mmread(os.path.join(hto_mat_folder, "matrix.mtx.gz" if args.gz == True else "matrix.mtx"))
features_path = os.path.join(hto_mat_folder, "features.tsv.gz" if args.gz == True else "features.tsv")
feature_ids = [row[0] for row in csv.reader(open(features_path), delimiter="\t")]
gene_names = [row[1] for row in csv.reader(open(features_path), delimiter="\t")]
barcodes_path = os.path.join(hto_mat_folder, "barcodes.tsv.gz" if args.gz == True else "barcodes.tsv")
barcodes = [row[0] for row in csv.reader(open(barcodes_path), delimiter="\t")]

if __name__ == '__main__':
    hto_mat = pd.DataFrame.sparse.from_spmatrix(hto_mat_mtx)
    hto_mat.columns = barcodes
    hto_mat.insert(loc=0, column="feature_id", value=feature_ids)
    hto_mat.insert(loc=1, column="gene_name", value=gene_names)
    
    hto_mat.to_csv("matrix_filtered.csv", index=False)


'''

parser.add_argument('--n_threads', help='Number of threads to use.',type=int, default=1)
parser.add_argument('--genome', help='Reference genome name. If not provided, we will infer it from the expression matrix file.')
parser.add_argument('--alpha', help='The Dirichlet prior concentration parameter (alpha) on samples.', type=float, default=0.0)
parser.add_argument('--min_Num_Genes', help='The Dirichlet prior concenration parameter on the background noise',type=float, default=1.0)
parser.add_argument('--min_Num_Umis', help='The Dirichlet prior concenration parameter on the background noise',type=float, default=1.0)
parser.add_argument('--min_signal_hashtag', help='Any cell/nucleus with less than min_signal hashtags from the signal will be marked as Negative.',type=float, default=10.0)
parser.add_argument('--randomState ', help='The random seed used in the KMeans algorithm to separate empty ADT droplets from others. ',type=int, default=0)
parser.add_argument('--generateDiagnosticPlots', help='Any cell/nucleus with less than min_signal hashtags from the signal will be marked as Negative.',type=float, default=10.0)
parser.add_argument('--output',  help='Output name.  All outputs will use it as the prefix.',default="demuxEM_res")


import csv
import gzip
import os
import scipy.io
 
# define MEX directory
matrix_dir = "/opt/sample345/outs/filtered_feature_bc_matrix"
# read in MEX format matrix as table
mat_filtered = scipy.io.mmread(os.path.join(matrix_dir, "matrix.mtx.gz"))
 
# list of transcript ids, e.g. 'ENSG00000187634'
features_path = os.path.join(matrix_dir, "features.tsv.gz")
feature_ids = [row[0] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]
 
# list of gene names, e.g. 'SAMD11'
gene_names = [row[1] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]
 
# list of feature_types, e.g. 'Gene Expression'
feature_types = [row[2] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]
 
# list of barcodes, e.g. 'AAACATACAAAACG-1'
barcodes_path = os.path.join(matrix_dir, "barcodes.tsv.gz")
barcodes = [row[0] for row in csv.reader(gzip.open(barcodes_path, mode="rt"), delimiter="\t")]

import pegasusio as io
from demuxEM.tools import *
import pegasus as pg
import argparse
import multiprocessing
from multiprocessing import freeze_support

parser = argparse.ArgumentParser(description='Parser for DemuxEM - Demultiplexing ')
parser.add_argument('--rna_data',  help= 'Input raw RNA expression matrix in 10x hdf5 format.')
parser.add_argument('--hto_matrix', help= 'HTO (antibody tag) count matrix in mtx format.')
parser.add_argument('--alpha', help='The Dirichlet prior concentration parameter (alpha) on samples.', type=float, default=0.0)
parser.add_argument('--alpha_noise', help='The Dirichlet prior concenration parameter on the background noise',type=float, default=1.0)
parser.add_argument('--tol', help='Threshold used for the EM convergence',type=float, default=1e-6)
parser.add_argument('--n_threads', help='Number of threads to use.',type=int, default=1)
parser.add_argument('--min_signal', help='Any cell/nucleus with less than min_signal hashtags from the signal will be marked as Negative.',type=float, default=10.0)

parser.add_argument('--output',  help='Output name',default="demuxEm.csv")


args = parser.parse_args()

umi_mat = args.rna_data
hto_mat = args.hto_matrix
alpha_val = args.alpha
alpha_noise_val = args.alpha_noise
min_signal_val = args.min_signal
tol_val = args.tol
n_threads_val = args.n_threads


#The input is read with 2 libraries, only pegasus and pegasus io for the mtx
#Anyways, both matrices are transformed to multimodal
hto_data = io.read_input(hto_mat)
umi_data = pg.read_input(input_file=umi_mat)
if __name__ == '__main__':
    multiprocessing.freeze_support()
    pg.qc_metrics(umi_data)
    pg.identify_robust_genes(umi_data)
    pg.estimate_background_probs(hto_data)
    pg.demultiplex(umi_data, hto_data,alpha=alpha_val, alpha_noise=alpha_noise_val,min_signal=min_signal_val,tol=tol_val,n_threads=n_threads_val)
    umi_data.obs.to_csv(args.output)
    

'''