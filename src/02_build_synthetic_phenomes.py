import os
import pathlib
import pandas as pd
import json
from smart_open import open
import yaml

# ------------------------------------------------------------------------------- #
# PART 2 - create_synthetic_mouse_phenomes #
# ------------------------------------------------------------------------------- #
# Specify your Gretel API key

from getpass import getpass
import pandas as pd
from gretel_client import configure_session, ClientConfig

base_path = pathlib.Path(os.getcwd().replace("/synthetics", ""))
data_path = base_path / 'mice_data_set' / 'data'
data_path

pd.set_option('max_colwidth', None)

configure_session(ClientConfig(api_key=getpass(prompt="Enter Gretel API key"), 
                               endpoint="https://api.gretel.cloud"))
#

# model config file
with open( base_path / 'config' / "model_params.yml", 'r') as stream:
    config = yaml.safe_load(stream)

# print(json.dumps(config, indent=2))

# Define a function to submit a new model for a specific phenome batch dataset
def create_model(batch_num):
    seconds = int(time.time())
    project_name = "Training phenomes" + str(seconds)
    project = create_project(display_name=project_name)
    batchfile = "pheno_batch" + str(batch_num) + ".csv"
    trainpath = str(data_path / batchfile)
    model = project.create_model_obj(model_config=config)
    model.data_source = trainpath
    model.submit(upload_data_source=True)  
    return(model)
#
# Submit all the phenome batches to train in parallel; poll for completion

from gretel_client.helpers import poll
from gretel_client import create_project
import time

# Create a model for each batch
models = []
for i in range(7):
    model = create_model(i)
    models.append(model)

# Poll for completion. Resubmit errors.
training = True
while training:
    time.sleep(60)
    training = False
    print()
    for i in range(7):
        model = models[i]
        model._poll_job_endpoint()
        status = model.__dict__['_data']['model']['status']
        print("Batch " + str(i) + " has status: " + status)
        if ((status == "active") or (status == "pending")):
            training = True
        if status == "error":
            model = create_model(i)
            models[i] = model
            training = True           

# Now that models are complete, get each batches Synthetic Quality Score (SQS)            
batch = 0
print()
for model in models:
    model._poll_job_endpoint()
    status = model.__dict__['_data']['model']['status']
    if status == "error":
        print("Batch " + str(batch) + " ended with error")
    else:
        report = model.peek_report()
        sqs = report['synthetic_data_quality_score']['score'] # THIS GIVES ERROR WHEN NONE
        label = "Moderate"
        if sqs >= 80:
            label = "Excellent"
        elif sqs >= 60:
            label = "Good"
        print("Batch " + str(batch) + " completes with SQS: " + label + " (" + str(sqs) + ")")
    batch += 1

# Read the original phenome set

filename = "phenome_alldata.csv"
filepath = data_path / filename
phenome_orig = pd.read_csv(filepath)

phenome_orig.head()

# Merge the synthetic batches into one dataframe
# First gather the synthetic data for each batch

synth_batches = []
for i in range(7):
    model = models[i]
    synth = pd.read_csv(model.get_artifact_link("data_preview"), compression='gzip')
    synth_batches.append(synth)

# Merge batch 0 and 1 on common field sacweight
synth_batches[0]['g'] = synth_batches[0].groupby('sacweight').cumcount()
synth_batches[1]['g'] = synth_batches[1].groupby('sacweight').cumcount()
synth_allbatches = pd.merge(synth_batches[0],synth_batches[1],on=["sacweight", 'g'],how='left').drop('g', axis=1)

# Now merge in batch 2 on common fields SW16, SW20, SW17
synth_allbatches['g'] = synth_allbatches.groupby(['SW16','SW20', 'SW17']).cumcount()
synth_batches[2]['g'] = synth_batches[2].groupby(['SW16', 'SW20', 'SW17']).cumcount()
synth_allbatches = pd.merge(synth_allbatches,synth_batches[2],on=['SW16', 'SW20', 'SW17', 'g'],how='left').drop('g', axis=1)

# Now merge in batches 3 
synth_allbatches = pd.concat([synth_allbatches, synth_batches[3]], axis=1)

# Now merge in batch 4 using common fields 'methcage12', 'methcage9', 'methcage11', 'methcage8', 'methcage7', 'methcage10'
synth_allbatches['g'] = synth_allbatches.groupby(['methcage12', 'methcage9', 'methcage11', 'methcage8', 'methcage7', 'methcage10']).cumcount()
synth_batches[4]['g'] = synth_batches[4].groupby(['methcage12', 'methcage9', 'methcage11', 'methcage8', 'methcage7', 'methcage10']).cumcount()
synth_allbatches = pd.merge(synth_allbatches,synth_batches[4],on=['methcage12', 'methcage9', 'methcage11', 'methcage8', 'methcage7', 'methcage10', 'g'],how='left').drop('g', axis=1)

# Now merge in batch 5 using common fields 'methcage12', 'methcage9', 'methcage11', 'methcage8', 'methcage7', 'methcage10'
synth_allbatches['g'] = synth_allbatches.groupby(['methcage12', 'methcage9', 'methcage11', 'methcage8', 'methcage7', 'methcage10']).cumcount()
synth_batches[5]['g'] = synth_batches[5].groupby(['methcage12', 'methcage9', 'methcage11', 'methcage8', 'methcage7', 'methcage10']).cumcount()
synth_allbatches = pd.merge(synth_allbatches,synth_batches[5],on=['methcage12', 'methcage9', 'methcage11', 'methcage8', 'methcage7', 'methcage10', 'g'],how='left').drop('g', axis=1)

# Now merge in batches 6
synth_allbatches = pd.concat([synth_allbatches, synth_batches[6]], axis=1)
synth_allbatches

# Add back in the "id" and "discard" fields, and save off complete synthetic data

id_col = []
discard_col = []
for i in range(len(synth_allbatches.index)):
    id_col.append(i)
    discard_col.append("no")
    
synth_allbatches["id"] = id_col
synth_allbatches["discard"] = discard_col
filepath = data_path / 'phenome_alldata_synth.csv'
synth_allbatches.to_csv(filepath, index=False, header=True)

# Optional cell if you have already created the synthetic phenomes and just need a new seed file
# to analyze a new pheno. Be sure to set the data path at the top of the notebook first

import pandas as pd
filepath = data_path / 'phenome_alldata_synth.csv'
synth_allbatches = pd.read_csv(filepath)

# Save off the phenotypes values you plan to analyze so we can later condition the genotype synthesis
# with these values

filename = data_path / 'pheno_and_covariates.csv'
pheno_analysis_df = pd.read_csv(filename)
pheno_seeds = list(pheno_analysis_df["pheno_and_cov"])

print(pheno_seeds)

seeds_df = synth_allbatches.filter(pheno_seeds)
# The seeding won't work if there are any NaN's in the seedfile
seeds_df = seeds_df.fillna(0)
seeds_df.head()

len(seeds_df)

# When you create the seeds df, you must make sure that any rounding or casting to int
# is replicated when creating genome training files.

seedfile = data_path / 'phenome_seeds.csv'
seeds_df.to_csv(seedfile, index=False, header=True)

seeds_df.head()

# Generate report that shows the statistical performance between the training and synthetic data
# Use the synthetic batch that includes abBMD

from smart_open import open
from IPython.core.display import display, HTML


# Change batch_num to any value between 0 and 6 to view performance report for other batches
batch_num = 0
display(HTML(data=open(models[0].get_artifact_link("report")).read(), metadata=dict(isolated=True)))

