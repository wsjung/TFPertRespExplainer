import argparse
import sys
from os import listdir
import os.path
import numpy as np
import pandas as pd

def parse_args(argv):
    parser = argparse.ArgumentParser(description='')
    parser.add_argument(
        '-d', '--input_dir', required=True,
        help='Directory containing transgene induction files.')
    parser.add_argument(
        '-o', '--output_file', required=True,
        help='Output file.')
    parser.add_argument(
        '-t', '--hgnc_ensembl_table', required=True,
        help='HGNC to ENSEMBL table.') 
    # HGNC TO ENSEMBL CONVERSION TABLE FROM GENENAMES.ORG
    parsed = parser.parse_args(argv[1:])
    return parsed


# FULL GENE PD HEADER
###  tf_hgnc,tf_ensg,gene_ensg,log2FoldChange,padj

name_dict = {
    'tf': 'tf_hgnc',
    'ensembl': 'tf_ensg',
    'Gene': 'gene_ensg',
    'lfc': 'log2Foldchange',
    'padj': 'padj'
}


def main(argv):
    args = parse_args(argv)

    print('input_dir:', args.input_dir)
    print('output file:', args.output_file)
    print('HGNC to ENSEMBL table: ', args.hgnc_ensembl_table)

    tfs = sorted([f.split('.')[0] for f in listdir(args.input_dir) if os.path.isfile(os.path.join(args.input_dir, f))])

    # LOAD TABLE INTO PD
    table = pd.read_csv(args.hgnc_ensembl_table, sep='\t', header=0) # tab delimited
    table.columns = ['hgnc','ensemble']

    # CREATE EMPTY FULL GENE PD
    columns = ['tf_hgnc', 'tf_ensg', 'gene_ensg', 'log2FoldChange', 'padj']
    all_genes = pd.DataFrame(columns=columns)

    # LOOP AND LOAD EACH GENE INTO PD
    for tf in tfs:
        filename = tf + '.txt'
        print(filename)

        tf_df = pd.read_csv(os.path.join(args.input_dir, filename), sep='\t', header=0) # tab delimited

        # COUNT NUMBER OF ROWS IN GENE PD
        num_genes = len(tf_df.index)
        print(num_genes)

        # FIND CORRESPONDING ENSEMBL ID FROM TABLE
        ensembl = table.loc[table['hgnc'] == tf]['ensemble'].values[0]
        print(ensembl)

        # IF NO HGNC->ENSEMBL PAIR EXISTS --> PRINT AND SKIP
        if ensembl != ensembl:
            print('MISSING ENSEMBL ID FOR TF: ', tf)

        # ADD COLUMN OF ENSEMBL ID TO GENE PD
        tf_df['tf_ensg'] = ensembl

        print(tf_df.head())

        # COMPUTE PADJ AND ADD TO DATAFRAME
        ############ ASK SANDEEP, MICHAEL, AND YIMING: HOW IS P-ADJ CALCULATED? NUM_GENES * PVAL AND CLIPPING AT 1?
        tf_df['padj'] = tf_df['p_val'] * num_genes

        print(tf_df.head())

        # MAX OUT PADJ TO 1
        tf_df['padj'] = tf_df['padj'].clip(upper=1.0)

        print(tf_df.head())


        # ADD TF NAME AS COL
        tf_df['tf'] = tf

        print(tf_df.head())


        # DROP PVAL COL

        tf_df  = tf_df.drop(columns=['p_val'])

        print(tf_df.head())

        return

        # RENAME COLUMNS

        

        # APPEND GENE PD TO FULL GENE PD


        # EXPORT FULL GENE PD 

    with open(args.hgnc_ensembl_table, 'r') as table, open(args.output_file, 'w') as out:
        table.readline() # header

        line = table.readline()[:-1]
        while line:
            hgnc, ensembl = line.split('\t')

            if hgnc in files:
                if ensembl == '\n': # no ensembl id
                    continue

                out.write("%s\t%s\n" % (hgnc, ensembl))
                files.remove(hgnc) # remove converted genes

            line = table.readline()[:-1]

    print("missing ensemble ids: ", files)

if __name__ == '__main__':
    main(sys.argv)
