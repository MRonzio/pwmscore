# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 15:29:59 2024

@author: MRonzio

CCAAT-box score calculus
"""

import math
import pandas as pd
import argparse

# running example:
# CalcScoreJaspar.py -M MAA060_1_JASPAR2016.txt -o CTCAGCCAATCAGCGC -s pos


# COMMAND LINE OPTIONS
def mtoptions():
    parser = argparse.ArgumentParser(description='JASPAR matrix score calculator')    
    parser.add_argument('-M', '--matrix', dest='jaspar_matrix', required=True, help='jaspar PWM in jaspar format')
    parser.add_argument('-o', '--oligo', dest='test_oligo', required=True, help='target name')
    parser.add_argument('-s', '--strand', dest='strand', default="pos",
                         help='oligo strand. Options: "pos" for positive, "neg" for negative. Default "pos".')
    return parser.parse_args()

# jaspar_matrix="/home/ronziom/Scaricati/HumanMouse_ChIPseq22/data/double_CCAAT/MAA060_1_JASPAR2016"
# test_oligo="ACGACCCAATCAGCAG" # (test)
# ref_oligo_pos="CTCAGCCAATCAGCGC" # (best ma0060.1 for control)
# ref_oligo_neg="GCGCTGATTGGCTGAG" # (best ma0060.1 for control)


# Prepare matrix
def prep_mat(jaspar_matrix):
    # read PWM matrix files
    JASPAR_df = pd.read_csv(jaspar_matrix, header=None, sep=' ', index_col=0)

    # normalize jaspar
    pseudocount = 0.01
    matrix = (JASPAR_df/116 + pseudocount) / ((pseudocount)*4+1)
    return matrix

# Create negative matrix
def create_neg_matrix(mat):
    # reverse columns
    columns = mat.columns.tolist()
    columns = columns[::-1]
    neg_mat = mat[columns]
    # inverto valori delle righe, cambiando gli indici e poi ripristinando l'ordine
    neg_mat.rename(index={'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}, inplace=True)
    return neg_mat

# for each oligo this is the formula to calculate the score:
# score = (MSoligo - minM)/(maxM - minM)
# all factors are log transformed

# calculat max and min possible values for the matrix (i.e. maxM e minM)
def calcolo_max_min(mat, tipo):
    MS_score_list = []
    seq = ""
    product = 1
    for column in mat:
        if tipo == "max":
            nt = mat[column].idxmax()
        elif tipo == "min":
            nt = mat[column].idxmin()
        score_MS = mat.loc[nt, column]
        seq += nt
        MS_score_list.append(score_MS)
    for msscore in MS_score_list:
        product *= msscore
    return math.log(product)


def calcolo_score_da_oligo(mat, oligo, minM, limits):
    MS_score_list = []
    product = 1
    for column, letter in zip(mat, oligo):
        score_MS = mat.loc[letter, column]
        MS_score_list.append(score_MS)
    for msscore in MS_score_list:
        product *= msscore
    score = (math.log(product) - minM) / limits
    return score


def main():
    # parse options
    options = mtoptions()
    #print(f"Site : {set_options.site_ref}")
    site_ref=options.test_oligo
    strand=options.strand
    input_matrix=options.jaspar_matrix
    # read and prep matrix
    matrix=prep_mat(input_matrix)
    # calculate min and max for the choosen matrix
    maxM_value=calcolo_max_min(matrix, "max")
    minM_value=calcolo_max_min(matrix, "min")
    limits_value=maxM_value - minM_value
    oriented_matrix= create_neg_matrix(matrix) if strand=="neg" else matrix
    # calculate and output
    score=calcolo_score_da_oligo(oriented_matrix,oligo=site_ref,minM=minM_value,limits=limits_value)
    print(score)


if __name__ == '__main__':
    main() 
