import os
import pathlib
import pandas as pd
import json
from smart_open import open
import yaml

# ------------------------------------------------------------------------------- #
# PART 1 - build_phenome_training_data #
# ------------------------------------------------------------------------------- #
# SET DATA PATH
base_path = pathlib.Path(os.getcwd().replace("/synthetics", ""))
data_path = base_path / 'mice_data_set' / 'data'
data_path

# PARSE CONFIG FILE
with open( base_path / 'config' / "phenome_analysis.yml", 'r') as stream:
    analysis = yaml.safe_load(stream)
#
pheno_analysis = analysis['analysis'][0]['phenotypes']['params']['traits']

pheno_analysis_df = pd.DataFrame({"pheno": pheno_analysis})
filename = data_path / 'pheno_analysis.csv'
pheno_analysis_df.to_csv(filename, index=False, header=True)

pheno_posby_covariates = analysis['analysis'][0]['phenotypes']['params']['covariates']
pheno_and_cov_df = pd.DataFrame({"pheno_and_cov": pheno_posby_covariates})
filename = data_path / 'pheno_and_covariates.csv'
pheno_and_cov_df.to_csv(filename, index=False, header=True)

# MUSCLE AND BONE TRAITS AND COVARIATES
# ----------------------
# For all five muscle weights (TA, EDL, soleus, plantaris and
# gastrocnemius), we map QTLs conditioning on tibia length
# ("tibia"). For tibia length, we map QTLs conditioned on body weight.
#
# Tibia length explains 12-18% of variance in the muscle weights. The
# rationale for including tibia length as a covariate is bone length
# may somehow regulate muscle weight as well, and we would like to
# isolate the genetic factors that directly regulate development of
# the muscle tissues.
#  
# For bone-mineral density (BMD), we created a binary trait that
# signals "abnormal" BMD. We do not include any covariates when
# mapping QTLs for these traits. Note that body weight is also
# uncorrelated with BMD.
# 
# For all muscle and bone traits, we include a binary indicator for
# round SW16 as a covariate because the mice from this round showed
# substantial deviation in these traits compared to the rest of the
# mice.
bone_muscle_pheno = [
    'TA',
    'SW16',
    'tibia',
    'EDL',
    'soleus',
    'plantaris',
    'gastroc',
    'SW6',
    'sacweight',
    'BMD',
    'abBMD']

# OTHER PHYSIOLOGICAL TRAITS AND COVARIATES
# --------------------------
# Body weights bw1, bw2 and bw3 were measured on subsequent days of
# the methamphetamine sensitivity tests, and are highly correlated
# with each other (r^2 = 98%), so it is only necessary to map QTLs for
# one of them. The body weight measurements after sacrifice
# ("sacweight") show a considerable departure in Round SW17, so we
# include a binary indicator for this round as a covariate for
# sacweight. We include age as a covariate for the "bw0" body weight
# because it was measured while the mouse was still growing.
#
# Fasting glucose levels are explained partially by body weight (PVE =
# 6%), so we include body weight as a covariate for fasting glucose
# levels. Rounds SW1 and SW11 showed a considerable departure in
# fasting glucose levels from the other rounds, so we included binary
# indicators for these two rounds as covariates for fasting glucose
# levels.
other_physio_traits_pheno = [
    'bw0',
    'glucoseage',
    'bw1',
    'methage',
    'SW17',
    'PPIweight',
    'sacweight',
    'fastglucose',
    'SW1',
    'SW11',
    'taillength',
    'SW3',
    'SW4',
    'SW19',
    'SW20',
    'SW22',
    'SW24',
    'testisweight']

# FEAR CONDITIONING TRAITS AND COVARIATES
# ------------------------
# For all fear conditioning traits, the cage used for testing appears
# to have an effect on the phenotype, so we include binary indicators
# for cage as covariates for all FC phenotypes. Further, the FC
# phenotype measurements in Round SW17 show a noticeably different
# distribution in the FC phenotypes from the other rounds, so we
# include a binary indicator for round SW17 as a covariate in all FC
# traits.
#
# These analyses control for proportion of freezing on day 1 during
# exposure to the tone ("AvToneD1"). AvToneD1 explains 11-25% of the
# variance in the Day 2 and Day 3 freezing measures. Note that here we
# can map QTLs for freezing to the altered context on Day 3
# ("AvAltContextD3") as a quantitative trait after conditioning on
# AvToneD1 because the distribution for this trait is no longer quite
# so bimodal, and looks fairly "normal". So there is no need to map
# QTLs for the binary version of this trait.
#
# PreTrainD1 is a very ugly trait with massive box effects and a lot
# of low values, which might have to be removed as outliers. It is
# quite likely that these outliers represent the "deaf" mice that
# might be skewing the whole results. These outliers are present in
# every box, so not a box-specific effect.
fear_cond_traits_pheno = [
    'AVContextD2',
    'AVToneD1',
    'FCbox1',
    'FCbox2',
    'FCbox3',
    'SW17',
    'AVAltContextD3',
    'AvToneD3',
    'PreTrainD1',
    'SW10',
    'SW16',
    'SW20',
    'SW7',
    'SW14']

# METHAMPHETAMINE SENSITIVITY, LOCOMOTOR ACTIVITY AND ANXIETY-LIKE BEHAVIOR AND COVARIATES
# -------------------------------------------------------------------------
# We checked all the cages used in these tests to see whether the
# phenotypes measured using any given cage departed noticeably from
# the other cages. Cage #7 consistently has a large effect.
meth_loco_anxiety_pheno1 = [
    'D1totaldist0to15',
    'D1totaldist0to30',
    'D1TOTDIST5',
    'D1TOTDIST10',
    'D1TOTDIST15',
    'D1TOTDIST20',
    'D1TOTDIST25',
    'D1TOTDIST30',
    'D1ctrtime0to15',
    'D1ctrtime0to30',
    'D1hact0to15',
    'D1hact0to30',
    'D1vact0to15',
    'D1vact0to30',
    'methcage7',
    'methcage8',
    'methcage9',
    'methcage10',
    'methcage11',
    'methcage12']

meth_loco_anxiety_pheno2 = [
    'D2totaldist0to15',
    'D2totaldist0to30',
    'D2TOTDIST5',
    'D2TOTDIST10',
    'D2TOTDIST15',
    'D2TOTDIST20',
    'D2TOTDIST25',
    'D2TOTDIST30',
    'D2ctrtime0to15',
    'D2ctrtime0to30',
    'D2hact0to15',
    'D2hact0to30',
    'D2vact0to15',
    'D2vact0to30',
    'methcage7',
    'methcage8',
    'methcage9',
    'methcage10',
    'methcage11',
    'methcage12']

 

meth_loco_anxiety_pheno3 = [
    'D3totaldist0to15',
    'D3totaldist0to30',
    'D3TOTDIST5',
    'D3TOTDIST10',
    'D3TOTDIST15',
    'D3TOTDIST20',
    'D3TOTDIST25',
    'D3TOTDIST30',
    'D3ctrtime0to15',
    'D3ctrtime0to30',
    'D3hact0to15',
    'D3hact0to30',
    'D3vact0to15',
    'D3vact0to30',
    'methcage7',
    'methcage8',
    'methcage9',
    'methcage10',
    'methcage11',
    'methcage12']

# PREPULSE INHIBITION (PPI) PHENOTYPES AND COVARIATES
# ------------------------------------
# All boxes appear to have some effect on some of the PPI phenotypes,
# with Box #3 having a particularly large effect on some phenotypes,
# so we include all PPI box indicators as covariates in analysis of the
# PPI phenotypes.
#
# We also map QTLs for habituation to pulses by analyzing the startle
# response during the fourth block of pulse-alone trials against the
# startle response during the first block of pulse-alone trials.
ppi_pheno = [
    'pp3PPIavg',
    'pp6PPIavg',
    'pp12PPIavg',
    'PPIavg',
    'startle',
    'p120b4',
    'PPIbox1',
    'PPIbox2',
    'PPIbox3',
    'PPIbox4',
    'p120b1']


pheno_batches = [bone_muscle_pheno, other_physio_traits_pheno, fear_cond_traits_pheno, meth_loco_anxiety_pheno1, 
                 meth_loco_anxiety_pheno2, meth_loco_anxiety_pheno3, ppi_pheno]
#

# Read in the pheno data saved from map notebook and discard lines with no data

import pandas as pd
phenofile = data_path / "pheno_new.csv"
pheno_data = pd.read_csv(phenofile)
pheno_data = pheno_data[pheno_data["cageid"].notnull()]
pheno_data.head()

# Create synthetic training files for each phenome batch. Create one with "id" for joining with the genome
# data, and one without "id" for phenome training

for i in range(len(pheno_batches)):
    columns_use = pheno_batches[i]
    columns_use.append("id")
    pheno_batch = pheno_data.filter(columns_use) 
    pheno_batch.dropna(inplace=True)
    pheno_batch = pheno_batch.round(4)
    pheno_batch_file = "pheno_batch" + str(i) + "_withID.csv"
    filename = data_path / pheno_batch_file
    pheno_batch.to_csv(filename, header=True, index=False)
    pheno_batch = pheno_batch.drop(['id'], axis=1)
    pheno_batch_file = "pheno_batch" + str(i) + ".csv"
    filename = data_path / pheno_batch_file
    pheno_batch.to_csv(filename, header=True, index=False)
#
# Make one big phenome dataset using all the relevant columns

columns_use = []
for next_batch in pheno_batches:
    columns_use = columns_use + next_batch
    
# Remove duplicates
columns_uniq = list(set(columns_use))

# Filter down to just these columns
pheno_alldata = pheno_data.filter(columns_uniq)

# Save data out for later comparison with synthetic data
phenofile = data_path / "phenome_alldata.csv"
pheno_alldata.to_csv(phenofile, index=False, header=True)

