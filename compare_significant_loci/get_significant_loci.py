#!/usr/bin/env python3

import argparse
import gzip
import bgzip
import numpy as np
from intervaltree import Interval, IntervalTree
import pysam
import os

argparser = argparse.ArgumentParser(description = 'This script extract non-overlapping loci with significant hits from GWAS results.')


argparser.add_argument('-g', '--gwas', metavar = 'file', dest = 'in_gwas_file', type = str, required = True, help = 'Compressed (gzip) GWAS result file(s) in Regenie format.')
argparser.add_argument('-af', '--min-af', metavar = 'number', dest = 'min_af', type = float, required = False, default = 0, help = 'Threshold for the minimal alternate allele frequency. Default: 0')
argparser.add_argument('-q', '--min-imp-quality', metavar = 'float', dest = 'min_info', type = float, required = False, default = 0.0, help = 'Threshold for the minimal imputation quality (INFO or Rsq). Default: 0.0') 
argparser.add_argument('-o', '--output-dir', metavar = 'name', dest = 'output_dir', type = str, required = True, help = 'Phenotype name. Will be used as a prefix for the output files.')

#TODO : add allele frequency

CHROMOSOME_CODES = set(
        [str(i) for i in range(1, 24)] + 
        [f'chr{i}'for i in range(1, 23)] + 
        ['X', 'chrX']
    )

SIGNIFICANCE_THRESHOLD = 0.05 * 1e-6

def rename_chrom_hg38(chrom, pos):
    if chrom not in ['X', 'chrX', '23', '25']:
        return chrom
    if chrom in ['23', '25']:
        new_chrom = 'X'
    else:
        new_chrom = chrom
    if pos >= 10001 and pos <= 2781479:
        return f'{new_chrom}_par1'
    if pos >= 155701383 and pos <= 156030895:
        return f'{new_chrom}_par2'
    return f'{new_chrom}_nonpar'


def get_window_hg38(chrom, pos, window_size):
    if chrom.endswith('X_nonpar'):
        return max(pos - window_size, 2781480), min(pos + window_size, 155701382)
    elif chrom.endswith('X_par1'):
        return max(pos - window_size, 10001), min(pos + window_size, 2781479)
    elif chrom.endswith('X_par2'):
        return max(pos - window_size, 155701383), min(pos + window_size, 156030895)
    else:
        return max(pos - window_size, 0), pos + window_size


def read_regenie(gwas_filename, min_af, min_info, window_size):
    required_columns = ['chrom', 'pos', 'af', 'n_samples',  'pval']
    loci_by_chrom = dict()
    n_used = 0
    with gzip.open(gwas_filename, 'rt') as ifile:
        print(f'Scanning {gwas_filename}...')
        header = ifile.readline().rstrip().split('\t')
        if any(c not in header for  c in required_columns):
            raise Exception('Required chrom, pos, af, n_samples, and pval columns are missing from the GWAS header.')
        chrom_idx = header.index('chrom')
        genpos_idx = header.index('pos')
        a1freq_idx = header.index('af')
        n_idx = header.index('n_samples')
        log10p_idx = header.index('pval')
        if 'INFO' in header:
            info_idx = header.index('INFO')
        else:
            info_idx = None
        for n, line in enumerate(ifile, 1):
            fields = line.rstrip().split('\t')
            if n % 1000000 == 0:
                print(f'\r{n} records scanned...', end = '', flush = True)
            chrom = fields[chrom_idx]
            if chrom not in CHROMOSOME_CODES:
                continue
            a1freq = float(fields[a1freq_idx])
            n_samples = int(fields[n_idx])
            if a1freq < min_af: # will be a litte bit inflated when in non-par X, but we can accept this
                continue
            if info_idx is not None and float(fields[info_idx]) < min_info:
                continue
            n_used += 1
            log10p = float(fields[log10p_idx])
            if log10p < SIGNIFICANCE_THRESHOLD:
                pos = int(fields[genpos_idx])
                chrom = rename_chrom_hg38(chrom, pos)
                window_start, window_stop = get_window_hg38(chrom, pos, window_size)
                loci_by_chrom.setdefault(chrom, IntervalTree()).addi(window_start, window_stop)
        print(f'\rDone. {n} records scanned; {n_used} records used (AF >= {min_af}, INFO >= {min_info}).\t\t')
    print(f'Number of significant hits by chromosome:')
    for chrom, loci in loci_by_chrom.items():
        print(f' {chrom}: {len(loci)}')
    print(f'Number of non-overlapping significant loci by chromosome:')
    for chrom, loci in loci_by_chrom.items():
        loci.merge_overlaps()
        print(f' {chrom}: {len(loci)}')
    return loci_by_chrom


def write_loci(gwas_filename, min_af, min_info, loci_by_chrom, output_dir):
    records_by_loci = dict()
    n_extracted = 0

    os.makedirs(output_dir, exist_ok=True)

    with gzip.open(gwas_filename, 'rt') as ifile:
        print(f'Extracting loci from {gwas_filename}...')
        header = ifile.readline().rstrip().split('\t')
        chrom_idx = header.index('chrom')
        genpos_idx = header.index('pos')
        a1freq_idx = header.index('af')
        n_idx = header.index('n_samples')
        if 'INFO' in header:
            info_idx = header.index('INFO')
        else:
            info_idx = None
        for n, line in enumerate(ifile, 1):
            if n % 1000000 == 0:
                print(f'\r{n} records scanned...', end = '', flush = True)
            fields = line.rstrip().split('\t')
            chrom = fields[chrom_idx]
            pos = int(fields[genpos_idx])
            chrom = rename_chrom_hg38(chrom, pos)
            if chrom not in loci_by_chrom:
                continue
            a1freq = float(fields[a1freq_idx])
            n_samples = int(fields[n_idx])
            if a1freq < min_af:
                continue
            if info_idx is not None and float(fields[info_idx]) < min_info:
                continue
            loci = loci_by_chrom[chrom]
            for interval in loci.at(pos):
                n_extracted += 1
                interval_name = f'{chrom}_{interval.begin}_{interval.end}'
                records_by_loci.setdefault(interval_name, list()).append(dict(zip(header, fields)))
        print(f'\rDone. {n} records scanned. {n_extracted} records extracted.\t\t') 
    with open(os.path.join(output_dir, f'{gwas_filename.split("/")[-1]}_loci_list.txt'), 'wt') as ofile_loci_list:
        for locus, records in records_by_loci.items():
            locus_filename = os.path.join(output_dir, f'{gwas_filename.split("/")[-1]}_{locus}.gz')
            records = sorted(records, key = lambda r: int(r['pos']))
            with open(locus_filename, 'wb') as raw:
                with bgzip.BGZipWriter(raw) as ofile_locus:
                    for i, record in enumerate(records, 0):
                        if i == 0:
                            ofile_locus.write('\t'.join(header).encode('utf8'))
                            ofile_locus.write('\n'.encode('utf8'))
                        ofile_locus.write('\t'.join(record[h] for h in header).encode('utf8'))
                        ofile_locus.write('\n'.encode('utf8'))
            pysam.tabix_index(locus_filename, line_skip = 1, seq_col = 0, start_col = 1, end_col = 1, force = True)
            locus_chrom, locus_start, locus_end = locus.rsplit('_', 2)
            ofile_loci_list.write(f'{gwas_filename}\t{locus_chrom}\t{locus_start}\t{locus_end}\t{os.path.abspath(locus_filename)}\n')
            print(f'Writing {locus_filename}')
    print('Done.')
    

if __name__ == '__main__':
    args = argparser.parse_args()

    loci_by_chrom = read_regenie(args.in_gwas_file, args.min_af, args.min_info, 500000)
    #print(f"{loci_by_chrom=}")
    write_loci(args.in_gwas_file, args.min_af, args.min_info, loci_by_chrom, args.output_dir)