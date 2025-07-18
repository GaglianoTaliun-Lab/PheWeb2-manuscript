#!/usr/bin/env python3

# This script compares the significant loci between male and female GWAS results
# It reads the significant loci files and finds the overlapping loci between the two files

import pandas as pd
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
import argparse

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--male", type=str, required=True, help="Path to the male significant loci file")
    parser.add_argument("-f", "--female", type=str, required=True, help="Path to the female significant loci file")
    parser.add_argument("-c", "--combined", type=str, required=True, help="Path to the female significant loci file")
    args = parser.parse_args()

    # Read the significant loci files
    male_loci = read_sig_loci(args.male)
    female_loci = read_sig_loci(args.female)
    combined_loci = read_sig_loci(args.combined)
    print(f"{male_loci.shape=}")
    print(f"{female_loci.shape=}")
    print(f"{combined_loci.shape=}")

    overlap_df = loci_overlap2(male_loci, female_loci, 'male', 'female')
    print(f"{overlap_df.shape=}")
    
    group_counts = overlap_df['SET_LOCI_GROUP_ID'].value_counts()
    shared_groups = group_counts[group_counts == 2].index
    shared_df = overlap_df[overlap_df['SET'] == 'male,female']
    print(f"{shared_df.shape=}")
    overlap_df.to_csv("overlap_df.csv", index=False)
    venn_for_2(overlap_df)

    overlap_df3 = loci_overlap3(male_loci, female_loci, male_loci, 'male', 'female', 'combined')
    print(f"{overlap_df3.shape=}")
    overlap_df3.to_csv("overlap_df3.csv", index=False)
    #venn_for_3(overlap_df3)

def read_sig_loci(file_path: str) -> pd.DataFrame:
    '''
    Read the significant loci file and return a pandas DataFrame with the following columns:
    - 0: phenotype name
    - 1: chromosome
    - 2: start position
    - 3: end position
    '''
    df = pd.read_csv(file_path, sep="\t", usecols=[0,1,2,3], header=None).rename(columns={0: 'pheno', 1: 'chr', 2: 'start', 3: 'end'})
    df = df[~df['pheno'].isin(['AGE_NMBR_COM', 'ADM_DCS_AGE_COM', 'WHO_MPAG_AG_COM'])].reset_index(drop=True)
    df['ANALYSIS'] = 'PRIMARY'
    df = df.rename(columns={
        'pheno': 'PHENOTYPE',
        'chr': 'CHROM',
        'start': 'LOCUS_START',
        'end': 'LOCUS_STOP'
    })
    df['CHROM'] = df['CHROM'].replace('X_nonpar', 23)
    df['CHROM'] = df['CHROM'].replace('X_par1', 24)
    df['CHROM'] = df['CHROM'].replace('X_par2', 25)
    df['CHROM'] = df['CHROM'].astype(int)
    df['ID'] = df.apply(lambda row: f"{row['PHENOTYPE']}${row['CHROM']}${row['LOCUS_START']}${row['LOCUS_STOP']}", axis=1)

    return df


def get_df_id(*args: pd.DataFrame) -> list[str]:
    '''
    Get the unique identifier of rows for each DataFrame
    '''
    id = []
    for df in args:
        id.extend([f"{df['pheno'].iloc[i]}${df['chr'].iloc[i]}${df['start'].iloc[i]}${df['end'].iloc[i]}" for i in range(len(df))])
    return id


def find_overlapping_loci(male_loci: pd.DataFrame, female_loci: pd.DataFrame, id_list: list[str]) -> pd.DataFrame:
    '''
    Find the overlapping loci between male and female significant loci
    overlap is defined as start1 <= end2 and end1 >= start2
    '''
    # overlap_df = pd.DataFrame(columns=['pheno', 'chr', 'start_male', 'end_male', 'start_female', 'end_female', 'group'])
    # overlap_df['group'] = []
    overlap_rows = [] 

    for id in id_list:
        pheno = id.split('$')[0]
        chr = id.split('$')[1]
        start = int(id.split('$')[2])
        end = int(id.split('$')[3])
        
        male = male_loci[(male_loci['pheno'] == pheno) & (male_loci['chr'] == chr) & (male_loci['start'] <= end) & (male_loci['end'] >= start)].copy()
        female = female_loci[(female_loci['pheno'] == pheno) & (female_loci['chr'] == chr) & (female_loci['start'] <= end) & (female_loci['end'] >= start)].copy()

        group = []
        if not male.empty:
            group.append('male')
        if not female.empty:
            group.append('female')
        
        if len(group) > 1:
            group_str = 'both'
        else:
            group_str = group[0]

        if group:
            overlap_rows.append({
                'pheno': pheno,
                'chr': chr,
                'start': start,
                'end': end,
                'group': group_str
            })


    return pd.DataFrame(overlap_rows)

def loci_overlap2(df1_hits: pd.DataFrame, df2_hits: pd.DataFrame, label1: str, label2: str):
    '''
    Find the overlapping loci between two DataFrames
    did the same job as find_overlapping_loci, more rigorous
    '''
    remaining_loci = set()
    loci_groups = dict()
    df_overlap = set()

    # helper functions    
    def get_locus_key(row):
        return f'{row.PHENOTYPE}:{row.CHROM}:{row.LOCUS_START}:{row.LOCUS_STOP}'

    def get_loci_group_id(*dfs):
        if len(loci_groups) == 0:
            return 1
        keys = set()
        for df in dfs:
            keys.update(df.apply(lambda row: get_locus_key(row), axis=1))
        group_ids = set()
        for key in keys:
            if key in loci_groups:
                group_ids.add(loci_groups[key])
        assert len(group_ids) <= 1
        return max(loci_groups.values()) + 1 if len(group_ids) == 0 else group_ids.pop()
        
    def add_loci_group_id(group_id, *dfs):
        keys = set()
        for df in dfs:
            keys.update(df.apply(lambda row: get_locus_key(row), axis=1))
        for key in keys:
            loci_groups[key] = group_id
        
    def add_loci_to_overlap(df, set_name, loci_group_id, label):
        for _, row in df.iterrows():
            key = get_locus_key(row)
            if key in remaining_loci:
                remaining_loci.remove(key)
            df_overlap.add((
                row.ID, row.CHROM, row.LOCUS_START, row.LOCUS_STOP,
                row.PHENOTYPE, set_name, loci_group_id, label
            ))

    # collect all unique loci from PRIMARY analysis
    remaining_loci.update(df1_hits[df1_hits.ANALYSIS == 'PRIMARY'].apply(get_locus_key, axis=1))
    remaining_loci.update(df2_hits[df2_hits.ANALYSIS == 'PRIMARY'].apply(get_locus_key, axis=1))

    # main loop
    while remaining_loci:
        phenotype, chrom, start, stop = remaining_loci.pop().split(':')
        chrom = int(chrom)
        start = int(start)
        stop = int(stop)

        df1_overlap = df1_hits[
            (df1_hits.ANALYSIS == 'PRIMARY') &
            (df1_hits.PHENOTYPE == phenotype) &
            (df1_hits.CHROM == chrom) &
            (df1_hits.LOCUS_START <= stop) &
            (df1_hits.LOCUS_STOP >= start)
        ]
        df2_overlap = df2_hits[
            (df2_hits.ANALYSIS == 'PRIMARY') &
            (df2_hits.PHENOTYPE == phenotype) &
            (df2_hits.CHROM == chrom) &
            (df2_hits.LOCUS_START <= stop) &
            (df2_hits.LOCUS_STOP >= start)
        ]

        if len(df1_overlap) >= 1 and len(df2_overlap) >= 1:
            loci_group_id = get_loci_group_id(df1_overlap, df2_overlap)
            add_loci_group_id(loci_group_id, df1_overlap, df2_overlap)
            add_loci_to_overlap(df1_overlap, f'{label1},{label2}', loci_group_id, label1)
            add_loci_to_overlap(df2_overlap, f'{label1},{label2}', loci_group_id, label2)
        elif len(df1_overlap) >= 1:
            loci_group_id = get_loci_group_id(df1_overlap)
            add_loci_group_id(loci_group_id, df1_overlap)
            add_loci_to_overlap(df1_overlap, f'{label1}', loci_group_id, label1)
        elif len(df2_overlap) >= 1:
            loci_group_id = get_loci_group_id(df2_overlap)
            add_loci_group_id(loci_group_id, df2_overlap)
            add_loci_to_overlap(df2_overlap, f'{label2}', loci_group_id, label2)

    return pd.DataFrame(df_overlap, columns=[
        'ID', 'CHROM', 'LOCUS_START', 'LOCUS_STOP', 'PHENOTYPE',
        'SET', 'SET_LOCI_GROUP_ID', 'SET_LOCUS_SOURCE'
    ])

def loci_overlap3(df1_hits: pd.DataFrame, df2_hits: pd.DataFrame, df3_hits: pd.DataFrame, label1: str, label2: str, label3: str) -> pd.DataFrame:
    '''
    Find the overlapping loci between three DataFrames
    '''

    remaining_loci = set()
    loci_groups = dict()
    df_overlap = set()

    # helper functions    
    def get_locus_key(row):
        return f'{row.PHENOTYPE}:{row.CHROM}:{row.LOCUS_START}:{row.LOCUS_STOP}'

    def get_loci_group_id(*dfs):
        if len(loci_groups) == 0: return 1
        keys = set()
        for df in dfs: keys.update(df.apply(lambda row: get_locus_key(row), axis = 1))
        group_ids = set()
        for key in keys:
            if key in loci_groups: group_ids.add(loci_groups[key])
        assert len(group_ids) <= 1
        return max(loci_groups.values()) + 1 if len(group_ids) == 0 else group_ids.pop()
        
    def add_loci_group_id(group_id, *dfs):
        keys = set()
        for df in dfs: keys.update(df.apply(lambda row: get_locus_key(row), axis = 1))
        for key in keys: loci_groups[key] = group_id
        
    def add_loci_to_overlap(df, set_name, loci_group_id, label):
        for index, row in df.iterrows():
            key = get_locus_key(row)
            if key in remaining_loci: remaining_loci.remove(key) # do not need to process this locus in the future
            df_overlap.add((row.ID, row.CHROM, row.LOCUS_START, row.LOCUS_STOP, row.PHENOTYPE, set_name, loci_group_id, label))

    # get all unique loci:
    remaining_loci.update(df1_hits[(df1_hits.ANALYSIS == 'PRIMARY')].apply(lambda row: get_locus_key(row), axis = 1))
    remaining_loci.update(df2_hits[(df2_hits.ANALYSIS == 'PRIMARY')].apply(lambda row: get_locus_key(row), axis = 1))
    remaining_loci.update(df3_hits[(df3_hits.ANALYSIS == 'PRIMARY')].apply(lambda row: get_locus_key(row), axis = 1))
    
    # iterate over all loci and check for occurances in multiple data frames
    # for locus in all_loci:
    while len(remaining_loci) > 0:
        phenotype, chrom, start, stop = remaining_loci.pop().split(':')
        chrom = int(chrom) # chromosomes coded as 1-23 in Regenie output
        start = int(start)
        stop = int(stop)
        df1_overlap = df1_hits[(df1_hits.ANALYSIS == 'PRIMARY') & (df1_hits.PHENOTYPE == phenotype) & 
                        (df1_hits.CHROM == chrom) & (df1_hits.LOCUS_START <= stop) & (df1_hits.LOCUS_STOP >= start)]
        df2_overlap = df2_hits[(df2_hits.ANALYSIS == 'PRIMARY') & (df2_hits.PHENOTYPE == phenotype) & 
                        (df2_hits.CHROM == chrom) & (df2_hits.LOCUS_START <= stop) & (df2_hits.LOCUS_STOP >= start)]
        df3_overlap = df3_hits[(df3_hits.ANALYSIS == 'PRIMARY') & (df3_hits.PHENOTYPE == phenotype) & 
                        (df3_hits.CHROM == chrom) & (df3_hits.LOCUS_START <= stop) & (df3_hits.LOCUS_STOP >= start)]
        if (len(df1_overlap) >= 1) and (len(df2_overlap) >= 1) and (len(df3_overlap) >= 1):
            loci_group_id = get_loci_group_id(df1_overlap, df2_overlap, df3_overlap)
            add_loci_group_id(loci_group_id, df1_overlap, df2_overlap, df3_overlap)
            add_loci_to_overlap(df1_overlap, f'{label1},{label2},{label3}', loci_group_id, label1)
            add_loci_to_overlap(df2_overlap, f'{label1},{label2},{label3}', loci_group_id, label2)
            add_loci_to_overlap(df3_overlap, f'{label1},{label2},{label3}', loci_group_id, label3)
        elif (len(df1_overlap) >= 1) and (len(df2_overlap) >= 1) and (len(df3_overlap) == 0):
            loci_group_id = get_loci_group_id(df1_overlap, df2_overlap)
            add_loci_group_id(loci_group_id, df1_overlap, df2_overlap)
            add_loci_to_overlap(df1_overlap, f'{label1},{label2}', loci_group_id, label1)
            add_loci_to_overlap(df2_overlap, f'{label1},{label2}', loci_group_id, label2)
        elif (len(df1_overlap) >= 1) and (len(df2_overlap) == 0) and (len(df3_overlap) >= 1):
            loci_group_id = get_loci_group_id(df1_overlap, df3_overlap)
            add_loci_group_id(loci_group_id, df1_overlap, df3_overlap)
            add_loci_to_overlap(df1_overlap, f'{label1},{label3}', loci_group_id, label1)
            add_loci_to_overlap(df3_overlap, f'{label1},{label3}', loci_group_id, label3)
        elif (len(df1_overlap) == 0) and (len(df2_overlap) >= 1) and (len(df3_overlap) >= 1):
            loci_group_id = get_loci_group_id(df2_overlap, df3_overlap)
            add_loci_group_id(loci_group_id, df2_overlap, df3_overlap)
            add_loci_to_overlap(df2_overlap, f'{label2},{label3}', loci_group_id, label2)
            add_loci_to_overlap(df3_overlap, f'{label2},{label3}', loci_group_id, label3)
        elif (len(df1_overlap) >= 1) and (len(df2_overlap) == 0) and (len(df3_overlap) == 0):
            loci_group_id = get_loci_group_id(df1_overlap)
            add_loci_group_id(loci_group_id, df1_overlap)
            add_loci_to_overlap(df1_overlap, f'{label1}', loci_group_id, label1)
        elif (len(df1_overlap) == 0) and (len(df2_overlap) >= 1) and (len(df3_overlap) == 0):
            loci_group_id = get_loci_group_id(df2_overlap)
            add_loci_group_id(loci_group_id, df2_overlap)
            add_loci_to_overlap(df2_overlap, f'{label2}', loci_group_id, label2)
        elif (len(df1_overlap) == 0) and (len(df2_overlap) == 0) and (len(df3_overlap) >= 1):
            loci_group_id = get_loci_group_id(df3_overlap)
            add_loci_group_id(loci_group_id, df3_overlap)
            add_loci_to_overlap(df3_overlap, f'{label3}', loci_group_id, label3)
        else:
            raise Exception() # should never occur
    return pd.DataFrame(df_overlap, columns = ['ID', 'CHROM', 'LOCUS_START', 'LOCUS_STOP', 'PHENOTYPE', 'SET', 'SET_LOCI_GROUP_ID', 'SET_LOCUS_SOURCE'])

def venn_for_2(df: pd.DataFrame) -> None:
    '''
    Plot the venn diagram of the significant loci,
    save male-only and female-only loci, and count by phenotype
    '''
    only_male = set(df[df['SET'] == 'male']['SET_LOCI_GROUP_ID'])
    only_female = set(df[df['SET'] == 'female']['SET_LOCI_GROUP_ID'])
    both = set(df[df['SET'] == 'male,female']['SET_LOCI_GROUP_ID'])

    shared = both
    male_only = only_male - shared
    female_only = only_female - shared

    # Filter full rows from df
    df_male_only = df[df['SET_LOCI_GROUP_ID'].isin(male_only)]
    df_female_only = df[df['SET_LOCI_GROUP_ID'].isin(female_only)]

    # Save to CSV
    df_male_only.to_csv('male_only_rows.csv', index=False)
    df_female_only.to_csv('female_only_rows.csv', index=False)

    # Count by PHENOTYPE
    count_by_phenotype = pd.concat([
        df_male_only.assign(SPECIFICITY='Male'),
        df_female_only.assign(SPECIFICITY='Female')
    ]).groupby(['PHENOTYPE', 'SPECIFICITY']).size().reset_index(name='COUNT')

    # Save counts to CSV (optional)
    count_by_phenotype.to_csv('sex_specific_loci_counts_by_phenotype.csv', index=False)

    venn2(subsets=(len(male_only), len(female_only), len(shared)),
      set_labels=('Male', 'Female'))

    plt.title("Fig1: Venn Diagram of Shared and Sex-specific Loci")
    plt.savefig("sex_specific_loci_venn.png", dpi=300)
    plt.show()

# def venn_for_3(df: pd.DataFrame) -> None:
#     '''
#     Plot the venn diagram of the significant loci
#     '''
#     only_male = set(df[df['SET'].isin(['male,none', 'none,male'])]['SET_LOCI_GROUP_ID'])
#     only_female = set(df[df['SET'] == 'female']['SET_LOCI_GROUP_ID'])
#     both = set(df[df['SET'] == 'male,female,none']['SET_LOCI_GROUP_ID'])

#     shared = both
#     male_only = only_male - shared
#     female_only = only_female - shared
#     venn2(subsets=(len(male_only), len(female_only), len(shared)),
#       set_labels=('Male', 'Female', 'Combined'))

#     plt.title("Fig2: Venn Diagram of Shared and Sex-specific Loci with three inputs")
#     plt.savefig("sex_specific_loci_venn2.png", dpi=300)
#     plt.show()

if __name__ == "__main__":
    main()