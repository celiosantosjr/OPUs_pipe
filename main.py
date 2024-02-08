import os
import gzip
import argparse
import pandas as pd   
from glob import glob 
from lib.utils import create_directory, corruption_test
from lib.internal_engine import loadprots, clean_records
from lib.internal_engine import cluster, getannotations


def setargs():

    parser = argparse.ArgumentParser()
    
    parser.add_argument('-ofolder', type=str, default='res', help='Output folder path')
    parser.add_argument('-pfolder', type=str, required=True, help='Project folder path')
    parser.add_argument('-annofolder', type=str, help='Annotation folder path')
    parser.add_argument('-annoxt', type=str, default="_emapper.annotations.tsv.gz", help='Annotation file extension')
    parser.add_argument('-minseqid', type=float, default=0.97, help='Minimum sequence identity')
    parser.add_argument('-threads', type=int, default=3, help='Number of threads')
    parser.add_argument('-minocc', type=int, default=2, help='Minimum number of occurrences to keep a pOTU')
    parser.add_argument('-minlen', type=int, default=98, help='Minimum number of residues to consider true a protein')
    parser.add_argument('-xt', type=str, default='.faa.gz', help='Number of threads')
    parser.add_argument('-maxmem', type=str, default='10G', help='Maximum memory (10G default)')
    
    return parser

    
def main():
    
    parser = setargs()    
    args = parser.parse_args()
    
    minocc = args.minocc
    minseqid = args.minseqid
    threads = args.threads
    annoxt = args.annoxt
    minlen = args.minlen
    maxmem = args.maxmem
    xt = args.xt
    
    pfolder = args.pfolder
    
    ofolder = args.ofolder
    create_directory(ofolder)
    
    if args.annofolder:
        annofolder = args.annofolder
    else:
        annofolder = f'{pfolder}/'
   
    for infile in glob(f'{pfolder}/*{xt}'):
    
        print(f"Working on file {infile}")
    
        sample = infile.split('/')[-1].replace(f'{xt}', '')
        oname_table = f'{ofolder}/{sample}_rename.tsv'
        oname_fasta = f'{ofolder}/{sample}_clean.fa.gz'
        annofile = f'{annofolder}/{sample}{annoxt}'
    
        c1, c2 = corruption_test(infile), corruption_test(annofile)
    
        if c1:
            print(f"File {infile} is corrupted")
            pass
        elif c2:
            print(f'File {annofile} is corrupted')
            pass
        else:
            records = loadprots(infile)
            clean_records(records,
                          oname_table=oname_table,
                          oname_fasta=oname_fasta,
                          minsize=minlen)
    
    print("Merging files")
    
    for infile in glob(f'{ofolder}/*.gz'):
     
        records = loadprots(infile)
        
        if len(records) == 0:
        
            print(f"The file {infile} is corrupted")
            pass
        
        else:    
        
            table = pd.DataFrame(records,
                                 columns=['seqname',
                                          'originalname',
                                          'sequence'])
            table = table[['originalname', 'sequence']]
            
            if len(table) > 0:
                
                with gzip.open(f'{ofolder}/summed_prots.faa.gz',
                               'at',
                               encoding='utf-8') as ofile:
                    
                    for _, h, s in table.itertuples():
                        ofile.write(f'>{h}\n{s}\n')
        
        os.remove(f'{infile}')    
    
    print("Wow, now we are clustering")
    
    pOTU_table = cluster(f'{ofolder}/summed_prots.faa.gz',
                         ofolder,
                         minseqid,
                         minocc,
                         threads,
                         maxmem)

    getannotations(pOTU_table,
                   ofolder,
                   annofolder,
                   annoxt)


if __name__ == '__main__':
    main()
