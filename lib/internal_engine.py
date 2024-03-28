import os
import gzip
import shutil
import pandas as pd

from Bio import SeqIO
from glob import glob
from subprocess import run
from collections import Counter
from lib.utils import create_directory
from lib.utils import emapper_ofile_fmt, cr
from lib.utils import corruption_test, mark_bad
from lib.utils import highest_value_with_percentage

def loadprots(file: str):
    '''
    It loads the fasta file containing proteins as
    a pandas dataframe.
    '''
    if file.endswith("xz") or file.endswith("gz"):
        if corruption_test(file):
            return "" 
    if file.endswith('gz'):
        file = gzip.open(file, 'rt', encoding='utf-8')
    records = list()
    nmb = 0
    for record in SeqIO.parse(file, "fasta"):
        newname = (record.id).split('.')[0]
        newname = newname + '_' + str(nmb) 
        records.append((newname, record.id, str(record.seq)))
        nmb += 1
    return records


def clean_records(records: list, oname_table: str, oname_fasta: str, minsize=98):
    '''
    It dereplicates sequences by identical length and sequence and 
    fragments as well, it outputs a fasta file of new sequences with
    completely new loads the fasta file containing proteins as
    a pandas dataframe.
    '''
    table = pd.DataFrame(records,
                         columns=['seqname',
                                  'originalname',
                                  'sequence'])
    table['good_seq'] = table.sequence.apply(lambda x: mark_bad(x, minsize))
    table.loc[table.good_seq == True,
              'sequence'] = table.loc[table.good_seq == True,
                                      'sequence'].apply(lambda x: x.upper().split('X')[0])
    print('Exporting renaming files')
    if not oname_table.endswith('xz'): oname_table += '.xz'    
    table.drop('sequence', axis=1).to_csv(oname_table,
                                          sep='\t', 
                                          header=True,
                                          index=None)
    print('Cleaning sequences')
    table = table[table.good_seq == True]
    print('Exporting fasta files')
    if not oname_fasta.endswith('gz'):
        oname_fasta += '.gz'
    with gzip.open(oname_fasta, 'wt', encoding='utf-8') as ofile:
        for _, seqname, protein in table[['seqname', 'sequence']].itertuples():
            ofile.write(f'>{seqname}\n{protein}\n')


def run_clustering(infile, ofolder, minseqid, threads, maxmem):
    run(['mmseqs',
         'easy-linclust',
         '--split-memory-limit', str(maxmem),
         '--cov-mode', '1',
         '--min-seq-id', str(minseqid),
         '--seq-id-mode', '1',
         '-e', '1e-5',
         '--threads', str(threads),
         infile,
         f'{ofolder}/result',
         'tmp/'])
    
    shutil.rmtree('tmp/')     


def process_cluster_table(cluster_df, minocc):
    print('# Extract sample information from sequence names')
    cluster_df['sample'] = cluster_df.sequence.apply(lambda x: x.split('_')[0])
    print('# Aggregate cluster information')
    keep_otus = cluster_df.representative.value_counts()
    keep_otus = keep_otus[keep_otus >= minocc].index
    OPU_table = cluster_df.groupby(['representative', 'sample']).size().reset_index(name='number')
    OPU_table = OPU_table[OPU_table.representative.isin(keep_otus)]
    print('# Pivot table for OPU representation')
    OPU_table = OPU_table.pivot_table(index='representative',
                                      columns='sample',
                                      values='number').fillna(0).astype('int').reset_index()
    print('# Add OPU column')
    OPU_table['OPU'] = ['OPU'+str(x) for x in OPU_table.index]
    print('# Reorder columns')
    OPU_table = OPU_table[['OPU', 'representative', *OPU_table.columns[1:-1]]]
    return OPU_table
    

def clean_up_files(ofolder):
    for x in ['result_all_seqs.fasta', 'result_cluster.tsv', 'result_rep_seq.fasta']:
        os.remove(f'{ofolder}/{x}')


def cluster(infile, ofolder, minseqid=0.97, minocc=2, threads=3, maxmem='10G'):
    print('# Making OPU table')
    run_clustering(infile, ofolder, minseqid, threads, maxmem)

    # Read cluster results
    cluster_df = pd.read_table(f'{ofolder}/result_cluster.tsv',
                               names=['representative', 'sequence']).drop_duplicates()

    # Save OPUs cluster relationship
    cluster_df.to_csv(f'{ofolder}/OPUs_cluster_relationship.tsv.xz',
                      sep='\t', header=True, index=None)
  
    # Save OPU table
    OPU_table = process_cluster_table(cluster_df, minocc)
    OPU_table.to_csv(f'{ofolder}/OPU_table.tsv.xz',
                      sep='\t', header=True, index=None)

    # Compress representative sequence file
    cr(f'{ofolder}/result_rep_seq.fasta', f'{ofolder}/result_rep_seq.fasta.xz')

    clean_up_files(ofolder)

    return OPU_table

  
def getannotations(df, ofolder, annofolder, annoxt):
    print('# Getting functions')

    print('# Load and merge dataframes')
    df2 = pd.read_table(f'{ofolder}/OPUs_cluster_relationship.tsv.xz')
    df2 = df2.merge(df[['representative', 'OPU']],
                    on='representative').rename(columns={'sequence': 'seqname'})

    print('# Process each namefile')
    for idx, namefile in enumerate(glob(f'{ofolder}/*_rename.tsv.xz')):
        x = pd.read_table(namefile)
        x = x[x.good_seq].merge(df2, on='seqname')[['originalname', 'OPU']]
        x = x.reset_index(drop=True)
        sample = x.loc[0, 'originalname'].split('.')[0]

        # Load and merge annotation file
        annofile = annofolder + sample + annoxt
        annofile = pd.read_table(annofile,
                                 names=emapper_ofile_fmt)[['originalname', 'KEGG_ko']]
        x = x.merge(annofile,
                    on='originalname',
                    how='left')

        # Save annotated file
        x.to_csv(f'{ofolder}/annotation_result_{idx}.tsv.xz',
                 sep='\t',
                 header=True,
                 index=None)

    print('# Concatenate annotated files')
    W = pd.concat([pd.read_table(infile) for infile in glob(f'{ofolder}/annotation_result_*.tsv.xz')]).fillna('UNKNOWN')
    W.to_csv(f'{ofolder}/general_OPUs_annotation.tsv.xz',
             sep='\t',
             header=True,
             index=None)

    print('# Summarize annotations')
    W_grouped = W.groupby(['OPU']).apply(lambda x: highest_value_with_percentage(dict(Counter(x.KEGG_ko))))
    W_grouped.reset_index().rename({0: 'annotation'},
                                   axis=1).to_csv(f'{ofolder}/summarized_OPUs_annotation.tsv.xz',
                                                  sep='\t',
                                                  header=True,
                                                  index=None)

    print('# Remove temporary files')
    [os.remove(infile) for infile in glob(f'{ofolder}/annotation_result_*.tsv.xz')]
    
