import pandas as pd
import argparse
import tomllib
import numpy as np
from scipy.stats import skew, kurtosis, rankdata, norm
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
import os
import re

argparser = argparse.ArgumentParser(description = 'Prepares summary stats about continuous phenotypes using a given csv file and a dictionarry.')
argparser.add_argument('-c', '--config', metavar = 'name', dest = 'config_file', type = str, required = True, help = "Config file in TOML format containing file location, sheet name and relevant column names. \n     For more information, see README")
argparser.add_argument('-o', '--output', metavar = 'name', dest = 'out_prefix', type = str, required = True, help = 'Prefix for output files.')
argparser.add_argument('-t', '--outlier_threshold', metavar = 'value', dest = 'outlier_value', type = int, default=5, help = 'Standard threshold (number of fold above and below of threshold) Fold value of the given outlier detection threshold (default 10)')

args = argparser.parse_args()

with open(args.config_file, mode="rb") as fp:
    TOML = tomllib.load(fp)

COLS_CODE   = TOML['COLUMN_NAMES_CODES']
COLS_MAIN   = TOML['COLUMN_NAMES_MAIN']
COLS_PHENO  = TOML['PHENOTYPE_COLS']
MISSING     = TOML['MISSING']

#these filters use the dict_categories values, but it would be better if they used the unique values from the actual data
#remove all the variables we want to ignore, if there's only 2 left, it is a binary variable
def filter_binary_variables(df, variable_missing_codes):
    label_cols = df.columns.values.tolist()
    for label in label_cols:
        df_codes = df[label].values.tolist()
        if (variable_missing_codes.get(label) != None):
            df_only_codes = np.unique(list(map(lambda x: np.nan if x in variable_missing_codes.get(label) else x, df_codes)))
        else:
            df_only_codes = np.unique(df_codes)
        df_only_codes = [x for x in df_only_codes if x == x]
        if len(df_only_codes) != 2:
            continue
        yield label
    return

#remove all the variables we want to ignore, if there's more than 2 left, it is a categorical (or quantitative) variable
#correction, lets do less than 30 unique variables instead of 2 to get rid of semi quantitative.
def filter_categorical_variables(df, variable_missing_codes):
    label_cols = df.columns.values.tolist()
    for label in label_cols:
        df_codes = df[label].values.tolist()
        if (variable_missing_codes.get(label) != None):
            df_only_codes = np.unique(list(map(lambda x: np.nan if x in variable_missing_codes.get(label) else x, df_codes)))
        else:
            df_only_codes = np.unique(df_codes)
        df_only_codes = [x for x in df_only_codes if x == x]
        if len(df_only_codes) <= 2 or len(df_only_codes) > 20:
            continue
        count = len(df[label][df[label].isin(df_only_codes)])
        print(label, len(df_only_codes), count )
        yield label
        
    return

def populate_missing_codes(df):

    #get unique values in each variable name that is the same..
    dict = {'varname':[],
        'Missing codes':[]
       }

    missing = pd.DataFrame(dict)

    for name in df[COLS_CODE['VARIABLE']].unique():
        subset = df[df[COLS_CODE['VARIABLE']] == name]

        terms_to_avoid = '|'.join(MISSING['MISSING_CATEGORIES_SUB'])

        rows_missing_exact = subset[COLS_CODE['CATEGORY']].isin(MISSING["MISSING_CATEGORIES_EXACT"])
        
        if terms_to_avoid != "":
            rows_missing_sub = subset[COLS_CODE['CATEGORY']].str.contains(terms_to_avoid)
        else:
            rows_missing_sub = subset[COLS_CODE['CATEGORY']].str.contains(terms_to_avoid) == False

        missing_codes = subset[ rows_missing_exact | rows_missing_sub ][COLS_CODE["CODE"]].astype(str)

        codes = ",".join(missing_codes)
        missing.loc[len(missing.index)] = [name, codes] 
    return missing

def filter_continuous_variables(df_variables, variable_missing_codes, df_pheno, df_pheno_final):

    for index, row in df_variables.iterrows():

        print(f"starting {index} out of {len(df_variables)}", flush = True)

        variable = row[COLS_MAIN['VARNAME']]
        print(variable, flush = True)
        if variable not in df_pheno.columns:
            print("not in columns", flush = True)
            continue
        if variable in df_pheno_final.columns:
            print("already in the final (duplicate)", flush = True)
            continue
        if (df_pheno.dtypes[variable] != 'float64' and df_pheno.dtypes[variable] != 'int64'):
            print("needs to be a float or int", flush = True)
            continue
        if 'RX_MED' in variable:
            print("RX_MED not relevant for GWAS", flush = True)
            continue

        unit = row[COLS_MAIN['UNITTYPE']]
        missing_codes = variable_missing_codes.get(variable, set())

        df_pheno_cut = df_pheno[[COLS_PHENO['SEX_VAR'], variable]].copy()
        
        #Removes Known missing codes
        df_pheno_cut['RECODED'] = df_pheno_cut[variable].apply(lambda x: None if x in missing_codes else x )

        df_not_missing = df_pheno_cut[~df_pheno_cut.RECODED.isna()]
        if len(df_not_missing) == 0:
            continue

        mean = np.nanmean(df_pheno_cut.RECODED)
        unique_values = len(df_not_missing.RECODED.unique())
        ## Records several assessement of Outliers
        sd_details={}
        sd=np.nanstd(df_pheno_cut.RECODED)

        if sd != 0 :
            for i in range(2,int(args.outlier_value)+1):
                sd_details['[Outliers] Number of values below mean -'+str(i)+' SD'] = len(df_not_missing.RECODED[df_not_missing.RECODED < (mean-i*sd)])
                sd_details['[Outliers] Number of values above mean +'+str(i)+' SD'] = len(df_not_missing.RECODED[df_not_missing.RECODED > (mean+i*sd)])
                if sd_details['[Outliers] Number of values below mean -'+str(i)+' SD'] !=0 :
                    sd_details['[Outliers] Maximum value below mean -'+str(i)+' SD'] = max(df_not_missing.RECODED[df_not_missing.RECODED < (mean-i*sd)])
                else :
                    sd_details['[Outliers] Maximum value below mean -'+str(i)+' SD'] = 'NA'
                if sd_details['[Outliers] Number of values above mean +'+str(i)+' SD'] !=0 :
                    sd_details['[Outliers] Minimal value above mean +'+str(i)+' SD'] = min(df_not_missing.RECODED[df_not_missing.RECODED > (mean+i*sd)])
                else :
                    sd_details['[Outliers] Minimal value above mean +'+str(i)+' SD'] = 'NA'
            min_value = df_not_missing.RECODED.min()
            max_value = df_not_missing.RECODED.max()
            PERC_5,PERC95 = np.percentile(df_not_missing.RECODED , [5,95])
            S=skew(df_not_missing.RECODED, axis=0, bias=True)
        else :
            ## if SD=0 only one value is present therefor these value are irrelevant
            PERC_5='NA'
            PERC95='NA'
            S='NA'

        #If the males females are string (such as "M", "F" or "male", "female") we need to recode them to binary.
        # do this before!! not at each loop
        if (isinstance(COLS_PHENO['MALE'], str)): 
            df_not_missing.loc[:,COLS_PHENO['SEX_VAR']] = df_not_missing.loc[:,COLS_PHENO['SEX_VAR']].replace(COLS_PHENO['MALE'], 0)
        if (isinstance(COLS_PHENO['FEMALE'], str)): 
            df_not_missing.loc[:,COLS_PHENO['SEX_VAR']] = df_not_missing.loc[:,COLS_PHENO['SEX_VAR']].replace(COLS_PHENO['FEMALE'], 1)

        plot(df_not_missing, row, mean, sd) 

        # Several Statistics
        min_value = df_not_missing.RECODED.min()
        max_value = df_not_missing.RECODED.max()
        df_males = df_not_missing[df_not_missing[COLS_PHENO['SEX_VAR']] == 0]
        df_females = df_not_missing[df_not_missing[COLS_PHENO['SEX_VAR']] == 1]
        mean_value = df_not_missing.RECODED.mean()
        median_value = df_not_missing.RECODED.median()
        value_counts = df_not_missing.RECODED.value_counts()
        mode_freq = value_counts.values[0]
        mode_value = value_counts.index[0]
        n_total = len(df_not_missing)
        n_males = len(df_males)
        n_females = len(df_females)
        problem=[] ## Provides flags
        if n_total < 100 :
            problem.append('Total')
        if n_males == 0 :
            problem.append('Males')
        if n_females ==0 :
            problem.append('Females')
        if mode_freq / n_total > 0.10 : #if mode is more than 10% of values, flag problem.
            problem.append('[Statistics] Mode frequency')
        if unique_values < 30 :
            problem.append('[Statistics] N uniques values')

        print(f"done {index} out of {len(df_variables)}", flush = True)
        yield {**{
                '[Description] Variable' : variable, 
                '[Description] Domain':  row[COLS_MAIN['DATABASE']],
                '[Description] Label' : row[COLS_MAIN['LABEL']],
                '[Samples] Total': n_total,
                '[Samples] Males' : n_males,
                '[Samples] Females' : n_females,
                '[Statistics] Minimum': min_value,
                '[Statistics] Maximum' : max_value, 
                '[Statistics] N uniques values': unique_values,
                '[Statistics] PERC_5' : PERC_5,
                '[Statistics] median' : median_value,
                '[Statistics] PERC_95' : PERC95,
                '[Statistics] Mean': mean_value,
                '[Statistics] Mode': mode_value,
                '[Statistics] Mode frequency': mode_freq,
                '[Hidden] Skewness': S,
                '[Statistics] SD' : sd,
                '[Hidden] problem' : problem 
            },** sd_details,**{
                '[Description] Unit': unit,
                '[Description] Missing code' : '|'.join([str(i) for i in missing_codes])
            }}

def rint(listed): # Rank Inverse Normal Transformation
    c=3/8
    rank=rankdata(listed)
    x = (rank - c) / (len(rank) - 2*c + 1)
    return norm.ppf(x)

def plot(df,row,meaned,sd): # Creates three pannel displaying the phenotype distribution
    #Raw Phenotype Distribution
    maxVal=np.nanmax(df.RECODED)
    minVal=np.nanmin(df.RECODED)
    f, axes = plt.subplots(2, 3, height_ratios=[3,1],figsize=(33, 10))

    try:
        sns.histplot(data=df, ax=axes[0,0], x="RECODED")
    except np.core._exceptions._ArrayMemoryError:
        sns.kdeplot(data=df, ax=axes[0,0], x="RECODED")

    axes[0,0].axvline(x=meaned,color='black')

    if sd !=0 :
        ulim = (maxVal + sd) if (maxVal > args.outlier_value*sd+meaned) else (args.outlier_value*sd+sd+meaned)
        llim = (minVal - sd) if (minVal < meaned-args.outlier_value*sd) else (meaned-(args.outlier_value*sd)-sd)
    else :
        ulim = maxVal + 1
        llim = minVal - 1

    axes[0,0].set_xlim(llim, ulim)
    axes[1,0].boxplot(df.RECODED,vert=False)
    axes[1,0].set_xlim(llim, ulim)
    fold=2

    ## Displays SD outlier values
    if sd != 0 :
        while (fold*sd+meaned < ulim) :
            axes[0,0].axvline(x=fold*sd+meaned,color='red')
            fold+=1
            if fold >args.outlier_value:
                break
        fold=2
        while (meaned-(fold*sd) > llim) :
            axes[0,0].axvline(x=meaned-(fold*sd),color='red')
            fold+=1
            if fold > args.outlier_value:
                break

    axes[0,0].axes.get_xaxis().set_visible(False)
    axes[1,0].axes.get_yaxis().set_visible(False)
    plt.suptitle(f'Phenotype Distribution \n {row[COLS_MAIN["LABEL"]]}',size=17)
    axes[0,0].set_ylabel('Count')
    axes[1,0].set_xlabel(row[COLS_MAIN['UNITTYPE']])

    ## Normalized Values where appropriate
    if sd !=0 :
        normalized=rint(df.loc[:,'RECODED'])
        maxValn=np.nanmax(normalized)
        minValn=np.nanmin(normalized)
        try:
            sns.histplot(data=normalized, ax=axes[0,1])
        except np.core._exceptions._ArrayMemoryError:
            sns.kdeplot(data=normalized, ax=axes[0,1])
        axes[1,1].boxplot(normalized,vert=False)
        axes[1,1].set_xlim(minValn, maxValn)
        axes[0,1].set_xlim(minValn, maxValn)
        axes[0,1].set_ylabel('Count')
        axes[1,1].set_xlabel('Normalized Unit')
        axes[0,1].axes.get_xaxis().set_visible(False)
        axes[1,1].axes.get_yaxis().set_visible(False)
        separator=(maxVal-minVal)/100 # Display purposes
        separatorN=(maxValn-minValn)/100 # Display Purposes
        ## Assess Skewness
        S=skew(df.RECODED, axis=0, bias=True)
        SS=skew(df.loc[(df.loc[:,'RECODED']>meaned-(2*sd)) & (df.loc[:,'RECODED']<meaned+(2*sd)) ,'RECODED'], axis=0, bias=True)
        SSS=skew(df.loc[(df.loc[:,'RECODED']>meaned-(3*sd)) & (df.loc[:,'RECODED']<meaned+(3*sd)) ,'RECODED'], axis=0, bias=True)
        axes[0,0].annotate("Skewness = {S:.3f} \n Skewness (without outlier 2sd) = {SS:.3f}  \n Skewness (without outlier 3sd) = {SSS:.3f} \n ".format(
        S=S,
        SS=SS,
        SSS=SSS), xy=(minVal+separator, 1), xycoords='data', xytext=(0.01, .99), textcoords='axes fraction', va='top', ha='left', fontsize='medium')
    
    ## Provides QQ plot
    sm.qqplot(df.RECODED, line ='s',ax=axes[0,2])
    plt.suptitle(f'Phenotype Distribution \n {row[COLS_MAIN["LABEL"]]}',size=13)

    #create images directory to put the plots in.
    if not os.path.exists('images'):
        os.makedirs('images')

    #save figures to that directory.
    f.savefig(f'images/{row[COLS_MAIN["VARNAME"]]}_distribution.png', dpi=300)
    plt.close(f)

if __name__ == '__main__':

    #read in the excel spreadsheet
    df = pd.read_excel(TOML['DICT_FILE'], sheet_name = TOML['SHEET_NAME_CODES'])

    subset_values = []
    for value in COLS_CODE:
        if (COLS_CODE[value] != ""):
            subset_values.append(COLS_CODE[value])

    df = df[subset_values]
    df[COLS_CODE['CATEGORY']] = df[COLS_CODE['CATEGORY']].astype("string")

    #run function to create dataframe in needed format. E.g.:
    #   varname Missing codes
    # 0 ACUTE_RENAL_FAIL_ONSET_AGE  -9, -7, 77, 88, 99
    # 1 ADDICTION_DISORDER_ONSET_AGE    -9, -7, 77, 88, 99 
    df_missing_codes = populate_missing_codes(df)

    df_missing_codes = df_missing_codes[df_missing_codes.iloc[:,1] != ""]
    
    #remove all non digits, make all digits floats.
    variable_missing_codes = {row['varname']: set(None if re.search('[^0-9]+', x) else float(x) for x in row['Missing codes'].split(',')) for index, row in df_missing_codes.iterrows()}

    #read in phenotypes
    df_pheno = pd.read_csv(TOML['PHENOTYPE_FILE'], header = 0)

    #binary_variables = set(filter_binary_variables(df_pheno, variable_missing_codes)) 
    categorical_variables = set(filter_categorical_variables(df_pheno, variable_missing_codes))
    
    # categorical_variables)
    
    with open('categorical_variables', 'w') as f:
        for var in categorical_variables:
            f.write(f"{var}\n")
    

    # df_variables = pd.read_excel(TOML['DICT_FILE'], sheet_name = TOML['SHEET_NAME_MAIN'], header = 0)

    # df_variables[COLS_MAIN['UNITTYPE']] = df_variables[COLS_MAIN['UNITTYPE']].fillna('')

    # #remove binary and categorical variables (keeping only quantitative variables)
    # df_variables = df_variables[~df_variables[COLS_MAIN['VARNAME']].isin(binary_variables)]
    # df_variables = df_variables[~df_variables[COLS_MAIN['VARNAME']].isin(categorical_variables)]
    
    # #create final dataframe to help populate.
    # df_pheno = pd.read_csv(TOML['PHENOTYPE_FILE'], header = 0)
    # df_pheno_final = pd.DataFrame({'FID': df_pheno[COLS_PHENO['IID']], 'IID': df_pheno[COLS_PHENO['IID']]})

    # #need to make this more memory efficient
    
    # #populate variables for final dataframe.
    # print("Starting filter_continous variables...", flush = True)
    # df_variables_final = pd.DataFrame(filter_continuous_variables(df_variables, variable_missing_codes, df_pheno, df_pheno_final))

    # print("Done filtering continuous variables", flush = True)
    # df_variables_final.to_csv(f'{args.out_prefix}.summary.tsv', sep = '\t', header = True, index = False)
    # df_variables_final.to_json(f'{args.out_prefix}.summary.json', orient='records')
    # print("DONE!")