# ------------------------------------------------------------------------------- #
# PART 3 - build_genome_training_data #
# ------------------------------------------------------------------------------- #
import os
import pathlib
import pandas as pd

base_path = pathlib.Path(os.getcwd().replace("/synthetics", ""))
data_path = base_path / 'mice_data_set' / 'data' 
experiment_path = base_path / 'mice_data_set' / 'out' 
data_path

# Read in the geno data and remove the discards (slow)

import pandas as pd


genofile = data_path / "geno.txt"
geno = pd.read_csv(genofile, sep=' ')
geno = geno[geno["discard"] == "no"]
geno.head()

# Gather the subset of phenotype data you plan to analyze
# Note, it is important that any rounding you do on these values
# mirrors what you did when you created the seed file in 02

phenofile = data_path / "phenome_alldata.csv"

pheno_alldata = pd.read_csv(phenofile)

filename = data_path / 'pheno_analysis.csv'
pheno_analysis_df = pd.read_csv(filename)
pheno_analysis = list(pheno_analysis_df["pheno"])

filename = data_path / 'pheno_and_covariates.csv'
pheno_and_cov_df = pd.read_csv(filename)
pheno_and_cov = list(pheno_and_cov_df["pheno_and_cov"])

# Here we're going to set the SNP count to use per training batch to be
# 19 minus the number of phenotypes and covariates being analyzed
snp_cnt = 19 - len(pheno_and_cov)
print("Using " + str(snp_cnt) + " SNPs per training batch")

columns_use = pheno_and_cov
columns_use.append('id')
pheno = pheno_alldata.filter(columns_use).round(4)

# Now you must drop any NaN's in the subset of pheno's you decided to use
pheno = pheno.dropna()
pheno.head()


# Grab the original GWAS linear model results for each phenotype being analyzed
# These files also contain the chromosome position information

all_gwas_scores = []
for next_pheno in pheno_analysis:
    filename = "lm_" + next_pheno + "_1_79646.csv"
    gwasfile = experiment_path / filename
    gwas_scores = pd.read_csv(gwasfile)
    all_gwas_scores.append(gwas_scores)

all_gwas_scores[0].head()

# Sort each by pvalue

gwas_scores_sorted = []
for gwas_scores in all_gwas_scores:
    scores_sorted = gwas_scores.sort_values(by=['p']).reset_index()
    gwas_scores_sorted.append(scores_sorted)
    
gwas_scores_sorted[0].head()

# Group the SNPs into batches of snp_cnt (determined above)
# Set min_chromo and max_chromo if you plan to analyze all batches
# on certain chromozomes

interesting_threshold = 5e-8
#interesting_threshold = 5e-7
min_chromo = 11
max_chromo = 14

batches = {}
batch_num = -1
max_snps_per_batch = snp_cnt
snp_cnt = 0
last_chromo = -1
batch_avg_pvalue = {}
all_avg_pvalues = []
all_int_cnt = []
all_not_int_cnt = []
all_batch_nums = []
snp_to_batch = {}

# Get a list of the SNPs sort by position on the chromosome
chromo_sort = all_gwas_scores[0].sort_values(by=['chr','pos'], ascending=False).reset_index()

for i in range(len(chromo_sort)):
    chromo = chromo_sort.loc[i]['chr']
    pos = chromo_sort.loc[i]['pos'] 
    snp = chromo_sort.loc[i]['snp'] 
    if ((snp_cnt == max_snps_per_batch) or (chromo != last_chromo)):
        batch_num += 1
        batches[batch_num] = {}
        batches[batch_num]['chromos'] = []
        batches[batch_num]['chr_pos'] = []
        batches[batch_num]['snps'] = []
        snp_cnt = 0
        
    batches[batch_num]['chr_pos'].append(str(chromo) + "_" + str(pos))
    batches[batch_num]['chromos'].append(chromo)
    batches[batch_num]['snps'].append(snp)
    snp_to_batch[snp] = batch_num
    last_chromo = chromo
    snp_cnt += 1

# Now for each phenotype we're analyzing gather the pvalue interesting cnt per batch

pheno_batch_interesting = {}
for i, next_pheno in enumerate(pheno_analysis):
    pheno_batch_interesting[next_pheno] = {}
    for j in range(len(batches)):
        pheno_batch_interesting[next_pheno][j] = 0
        
    gwas_scores = gwas_scores_sorted[i]
    k = 0
    pscore = gwas_scores.loc[k]['p']
    snp = gwas_scores.loc[k]['snp']
    while pscore <= interesting_threshold:
        batch_num = snp_to_batch[snp]
        pheno_batch_interesting[next_pheno][batch_num] += 1
        k += 1
        pscore = gwas_scores.loc[k]['p']
        snp = gwas_scores.loc[k]['snp']
    
# Sum the interesting pvalues per batch

batch_interesting_sums = []
for i in range(len(batches)):
    interesting_sum = 0
    for next_pheno in pheno_analysis:
        interesting_sum += pheno_batch_interesting[next_pheno][i]
    batch_interesting_sums.append(interesting_sum) 

# Save the batches in the deemed interesting range

batches_in_chromo_range = []
for i in range(len(batches)):
    chromos = batches[i]['chromos']
    use_batch = False
    for chromo in chromos:
        if ((chromo >= min_chromo) & (chromo < max_chromo)):
            use_batch = True
    if use_batch:
        batches_in_chromo_range.append(i)
        

# Gather up the interesting and non-interesting batches

interesting_batches = [i for i in range(len(batches)) if batch_interesting_sums[i] > 0] 
noninteresting_batches = [i for i in range(len(batches)) if batch_interesting_sums[i] == 0]

# Optional cell to see which batches have interesting relationships to the pheno's you're studying

for i in range(len(batch_interesting_sums)):
    count = batch_interesting_sums[i]
    if count > 0:
        print("Batch " + str(i) + ": ")
        for next in pheno_analysis:
            cnt = pheno_batch_interesting[next][i]
            print("\t" + next + ": " + str(cnt))

# Create function to build a training set from a batch

def build_training(batch):
    
    training_min_rows = 25000
    
    # Gather the SNPs
    grp_columns = list(batches[batch]["snps"])
    grp_columns.append("id")
    geno_grp = geno.filter(grp_columns)

    # Round float values to integers
    floats = geno_grp.select_dtypes(include=['float64']) 
    for col in floats.columns.values:
        geno_grp[col] = round(geno_grp[col]).astype('int')    

    # Add in the phenome information to the genome training set
    genome_phenome = geno_grp.join(pheno.set_index('id'), on = "id", how = "inner")
    columns_use = list(geno_grp.columns)
    columns_use.reverse()
    columns_use = columns_use + pheno_and_cov  
    genome_train = genome_phenome.filter(columns_use)
    
    # Replicate training set to have a minimum of 25000 examples
    dataset_rows = len(genome_train)
    genome_train = pd.concat([genome_train] * (training_min_rows // dataset_rows + 1))
    genome_train.drop(['id'], axis=1, inplace=True)
    
    # Save the training file
    filename = "geno_batch" + str(batch) + "_train.csv"
    genofile = data_path / "genome_training_data" / filename
    genome_train.to_csv(genofile, index=False, header=True)
    
    # Now create of version of map.txt with just the SNPs in the training set of this first batch
    mapfile = data_path / "map.txt"
    mapdata = pd.read_csv(mapfile, sep=' ')
    mapdata_use = mapdata[mapdata["id"].isin(batches[batch]["snps"])]
    filename = "map_batch" + str(batch) + ".txt"
    mapfile_new = data_path / "genome_map_data" / filename
    mapdata_use.to_csv(mapfile_new, sep=' ', header=True, index=False)
    

# How many batches are in the intesesting chromosome range (if you chose to do that)
len(batches_in_chromo_range)

# If desired, create training sets for batches in the desired chromo range

batches_to_use = batches_in_chromo_range
all_batches_used = []
for batch in batches_to_use:
    build_training(batch)
    all_batches_used.append(batch)

# Alternatively, you can choose batches from the interesting and non-interesting batch lists.
# How many batches have interesting pvalues
len(interesting_batches)

# Create training sets for some or all of the batches with interesting pvalues
# Set the number you want to use below

import random

interesting_batches_to_use = 6
batches_to_use = random.sample(interesting_batches, interesting_batches_to_use)
all_batches_used = []
for batch in batches_to_use:
    build_training(batch)
    all_batches_used.append(batch)

# How many batches have interesting pvalues
len(noninteresting_batches)

# Also create training sets for some or all of the batches with no interesting pvalues

non_interesting_batches_to_use = 19911
batches_to_use = random.sample(noninteresting_batches, non_interesting_batches_to_use)
for batch in batches_to_use:
    build_training(batch)
    all_batches_used.append(batch)

# Save to a file the list of batch numbers we created training sets for

file_df = pd.DataFrame({"batch": all_batches_used})
filename = data_path / "batch_training_list.csv"
file_df.to_csv(filename, index=False, header=True)

