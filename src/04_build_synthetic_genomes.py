# ------------------------------------------------------------------------------- #
# PART 4 - create_synthetic_mouse_genomes #
# ------------------------------------------------------------------------------- #

# Specify your Gretel API key

from getpass import getpass
import pandas as pd
from gretel_client import configure_session, ClientConfig

pd.set_option('max_colwidth', None)

configure_session(ClientConfig(api_key=getpass(prompt="Enter Gretel API key"), 
                               endpoint="https://api.gretel.cloud"))

import os
import pathlib
import pandas as pd

base_path = pathlib.Path(os.getcwd().replace("/synthetics", ""))
data_path = base_path / 'mice_data_set' / 'data'
data_path

filename = data_path / "batch_training_list.csv"
training_list_df = pd.read_csv(filename)
batches = list(training_list_df["batch"])
batch_cnt = len(batches)
batch_cnt

import json
from smart_open import open
import yaml

# Read in the phenome seed df

seedfile = str(data_path / 'phenome_seeds.csv')
seed_df = pd.read_csv(seedfile)

# model config file
with open( base_path / 'config' / "model_params.yml", 'r') as stream:
    config = yaml.safe_load(stream)


# Add a seed task which will enable us to tie together the genome and phenome data

filename = data_path / 'pheno_analysis.csv'
pheno_analysis_df = pd.read_csv(filename)
pheno_analysis = list(pheno_analysis_df["pheno"])

filename = data_path / 'pheno_and_covariates.csv'
pheno_and_cov_df = pd.read_csv(filename)
seed_fields = list(pheno_and_cov_df["pheno_and_cov"])

task = {
    'type': 'seed',
    'attrs': {
        'fields': seed_fields
    }
}


# # Optimize parameters for complex dataset
# I WOULD TEST THE USER/DEFAULT FIRST
# config['models'][0]['synthetics']['task'] = task
# config['models'][0]['synthetics']['params']['epochs'] = 150
# config['models'][0]['synthetics']['params']['vocab_size'] = 0
# config['models'][0]['synthetics']['params']['rnn_units'] = 256
# config['models'][0]['synthetics']['params']['reset_states'] = False
# config['models'][0]['synthetics']['params']['learning_rate'] = 0.001
# config['models'][0]['synthetics']['privacy_filters']['similarity'] = None
# config['models'][0]['synthetics']['params']['dropout_rate'] = 0.5
# config['models'][0]['synthetics']['params']['gen_temp'] = 1.0
# config['models'][0]['synthetics']['generate']['num_records'] = len(seed_df)
# config['models'][0]['synthetics']['generate']['max_invalid'] = len(seed_df) * 10

seed_fields

seed_df.head()

len(seed_df)

from gretel_client import create_project

# Estabilish some variables that will enable the running of models in parallel

MAX_MODELS_RUNNING = 58 # I WOULD NEED TO INCLUDE THIS IN A CONFIG FILE
MAX_MODEL_CNT_IN_PROJECT = 8 # I WOULD NEED TO INCLUDE THIS IN A CONFIG FILE
total_model_cnt = batch_cnt
models_training = []
models_generating = []
models_complete = []
models_training_cnt = 0
models_generating_cnt = 0
models_complete_cnt = 0
models_error_cnt = 0
model_cnt_in_project = 0
project_to_model_mapping = {}
project_num = 0
batch_num = 0
base_project_name = "Illumina Genome Batch "
moreToDo = True
model_info = {}
gwas_errors = 0
gwas_error_batches = []

# Initialize first project

project_name = base_project_name + str(project_num)
current_project = create_project(display_name=project_name)


## Functions for training and generating multiple models in parallel

import time
def chk_model_completion():
    
    global models_training
    global models_generating
    global model_info
    global models_complete
    global models_training_cnt
    global models_generating_cnt
    global models_complete_cnt
    global models_error_cnt
    global model_cnt_in_project
    global project_num
    global current_project
    global project_to_model_mapping
    global project_name
    
    for model_name in models_training:
        
        model = model_info[model_name]["model"]
        try:
            model._poll_job_endpoint()
        except:
            model._poll_job_endpoint()
        status = model.__dict__['_data']['model']['status']
        print("Model " + model_name + " has training status: " + status)
        
        # If model ended in error or was lost, restart it
        if ((status == "error") or (status == "lost")):
            
            # Check if we need a new project
            if model_cnt_in_project >= MAX_MODEL_CNT_IN_PROJECT:
                project_num += 1
                project_name = base_project_name + str(project_num)
                try:
                    current_project = create_project(display_name=project_name)
                except:
                    current_project = create_project(display_name=project_name)
                model_cnt_in_project = 0
            
            # Start a new model
            
            filename = model_info[model_name]["filename"]
            filepath = data_path / "genome_training_data" / filename
            try:
                artifact_id = current_project.upload_artifact(filepath)
            except:
                artifact_id = current_project.upload_artifact(filepath)
            
            try:
                model = current_project.create_model_obj(model_config=config, data_source=artifact_id)
            except:
                model = current_project.create_model_obj(model_config=config, data_source=artifact_id)
            
            model_info[model_name]["model"] = model
            
            try:
                model.submit()
            except:
                model.submit()
                
            model_cnt_in_project += 1           
            print("Model restarted training due to error: " + str(model_name))
            models_error_cnt += 1
            if project_name in project_to_model_mapping:
                project_to_model_mapping[project_name]["models"].append(model_name)
            else:
                project_to_model_mapping[project_name] = {}
                project_to_model_mapping[project_name]["models"] = [model_name]
                project_to_model_mapping[project_name]["project"] = current_project
                project_to_model_mapping[project_name]["status"] = "active"
            
        # If completed, get SQS and start generating records from seeds
        if (status == 'completed'): 
            models_training.remove(model_name)
            models_training_cnt -= 1
            report = model.peek_report()
            if report:
                sqs = report['synthetic_data_quality_score']['score']
                model_info[model_name]["sqs"] = sqs  
                print("Model " + str(model_name) + " has SQS: " + str(sqs))
            
            # Generate more records using seeds
            
            print("Model started generating: " + model_name)
            try:
                rh = model.create_record_handler_obj(data_source=seedfile, params={"num_records": len(seed_df)})
                rh.submit_cloud()
            except:
                rh = model.create_record_handler_obj(data_source=seedfile, params={"num_records": len(seed_df)})
                rh.submit_cloud()
            model_info[model_name]["rh"] = rh
            models_generating.append(model_name)
            models_generating_cnt += 1
            
            

def chk_model_generation_completion():

    import numpy as np
    global models_generating
    global models_training
    global model_info
    global models_complete
    global models_generating_cnt
    global models_training_cnt
    global models_complete_cnt
    global models_error_cnt
    global model_cnt_in_project
    global project_num
    global current_project
    global project_to_model_mapping
    global project_name
    global gwas_errors
    global gwas_error_batches
    
    for model_name in models_generating:
        rh = model_info[model_name]["rh"]
        try:
            rh._poll_job_endpoint()
        except:
            rh._poll_job_endpoint()
        status = rh.__dict__['_data']['handler']['status']
        print("Model " + model_name + " has generating status: " + status)

        # If generation ends in error, restart by training a fresh model
        
        if ((status == "error") or (status == "lost")):
            
            # Check if we need a new project
            if model_cnt_in_project >= MAX_MODEL_CNT_IN_PROJECT:
                project_num += 1
                project_name = base_project_name + str(project_num)
                try:
                    current_project = create_project(display_name=project_name)
                except:
                    current_project = create_project(display_name=project_name)
                model_cnt_in_project = 0
            
            # Start a new model
            filename = model_info[model_name]["filename"]
            filepath = data_path / "genome_training_data" / filename
            try:
                artifact_id = current_project.upload_artifact(filepath)
            except:
                artifact_id = current_project.upload_artifact(filepath)
            
            
            try:
                model = current_project.create_model_obj(model_config=config, data_source=artifact_id)
            except:
                model = current_project.create_model_obj(model_config=config, data_source=artifact_id)

            model_info[model_name]["model"] = model
            try:
                model.submit()
            except:
                model.submit()
            model_cnt_in_project += 1           
            print("Model restarted training due to error in generation: " + str(model_name))
            models_error_cnt += 1
                  
            models_generating.remove(model_name)
            models_training.append(model_name)
            models_generating_cnt -= 1
            models_training_cnt += 1
            
            if project_name in project_to_model_mapping:
                project_to_model_mapping[project_name]["models"].append(model_name)
            else:
                project_to_model_mapping[project_name] = {}
                project_to_model_mapping[project_name]["models"] = [model_name]
                project_to_model_mapping[project_name]["project"] = current_project
                project_to_model_mapping[project_name]["status"] = "active"
            
        if status == "completed":
            
            models_generating.remove(model_name)
            models_generating_cnt -= 1
            models_complete.append(model_name)
            models_complete_cnt += 1
         
            synthetic_genomes = pd.read_csv(rh.get_artifact_link("data"), compression='gzip') 
            
            # Drop the phenome information from the genome synth data and add back in the fields "id" and "discard"

            id_col = []
            discard_col = []
            for i in range(len(synthetic_genomes.index)):
                id_col.append(i)
                discard_col.append("no")

            synthetic_genomes = synthetic_genomes.drop(seed_fields, axis=1)
            
            columns = ['id', 'discard']
            columns = columns + list(synthetic_genomes.columns)   
            synthetic_genomes["id"] = id_col
            synthetic_genomes["discard"] = discard_col
            synthetic_genomes = synthetic_genomes.filter(columns)
    
            # Save the synthetic data
        
            filename = "synthetic_genomes_" + model_name + ".txt"  
            filepath = data_path / "synthetic_genome_data" / filename
            synthetic_genomes.to_csv(filepath, index=False, sep=' ')  
            print("Synthetic data for " + model_name + " saved to: " + filename)
            
            # Compute an initial GWAS F1
            
            try:
                F1s = computeF1(model_name)
                print("GWAS F1 score fors " + model_name + " are: ")
                for i, next_pheno in enumerate(pheno_analysis):
                    print("\t" + next_pheno + ": " + str(F1s[i]))
                model_info[model_name]["F1"] = np.mean(F1s)
                    
            except:
                F1s = []
                print("GWAS error for model " + model_name)
                gwas_errors += 1
                gwas_error_batches.append(model_name)
                model_info[model_name]["F1"] = 0
                
            
            
                                                   
          

def chk_project_completion():
    
    global models_generating
    global models_training
    global project_to_model_mapping

    for next in project_to_model_mapping:
        if project_to_model_mapping[next]["status"] == "active":
            project_active = False
            for next_model in project_to_model_mapping[next]["models"]:
                if ((next_model in models_generating) or (next_model in models_training)):
                    project_active = True
            if (project_active == False):
                this_project = project_to_model_mapping[next]["project"]
                
                # Note the below line is what you'd comment out if you'd like to hold onto your
                # synthetic models and later run the synthetic performance report. This means you'll
                # have to manually go into the Gretel Console and delete unneeded projects, which can
                # grow quickly if you're processing all genomic batches
                
                this_project.delete()
                project_to_model_mapping[next]["status"] = "completed"
    

def start_more_models():
    
    import time
    
    global models_training
    global models_generating
    global model_info
    global models_training_cnt
    global models_generating_cnt
    global models_complete_cnt
    global model_cnt_in_project
    global current_project
    global filelist
    global batch_num
    global project_num
    global project_to_model_mapping
    global project_name
    
    while (((models_training_cnt + models_generating_cnt) < MAX_MODELS_RUNNING) and 
          ((models_training_cnt + models_generating_cnt + models_complete_cnt) < total_model_cnt)):
        
        # Check if we need a new project
        if model_cnt_in_project >= MAX_MODEL_CNT_IN_PROJECT:
            project_num += 1
            project_name = base_project_name + str(project_num)
            try:
                current_project = create_project(display_name=project_name)
            except:
                current_project = create_project(display_name=project_name)
            model_cnt_in_project = 0
            
        # Start a new model

        batch = batches[batch_num]
        batch_num += 1
        filename = "geno_batch" + str(batch) + "_train.csv"
        filepath = data_path / "genome_training_data" / filename
        df = pd.read_csv(filepath)
        cluster_size = len(df.columns)
        config['models'][0]['synthetics']['params']['field_cluster_size'] = cluster_size
        
        try:
            artifact_id = current_project.upload_artifact(filepath)
        except:
            artifact_id = current_project.upload_artifact(filepath)
          
        try:
            model = current_project.create_model_obj(model_config=config, data_source=artifact_id)
        except:
            model = current_project.create_model_obj(model_config=config, data_source=artifact_id)
        model_name = "batch" + str(batch)
        
        models_training.append(model_name)
        models_training_cnt += 1
        model_info[model_name] = {}
        model_info[model_name]["model"] = model
        model_info[model_name]["filename"] = filename
        try:
            model.submit()
        except:
            model.submit()
        model_cnt_in_project += 1
        print("Model started training: " + str(model_name))
        
        if project_name in project_to_model_mapping:
            project_to_model_mapping[project_name]["models"].append(model_name)
        else:
            project_to_model_mapping[project_name] = {}
            project_to_model_mapping[project_name]["models"] = [model_name]
            project_to_model_mapping[project_name]["project"] = current_project
            project_to_model_mapping[project_name]["status"] = "active"
            
        
        

def computeF1(model_name):
    
    from sklearn.metrics import f1_score

    base_path = pathlib.Path(os.getcwd().replace("/synthetics", ""))
    data_path = base_path / 'mice_data_set' / 'data' 
    real_gwas_path = base_path / 'mice_data_set' / 'out' 
    synthetic_gwas_path = base_path / 'mice_data_set' / 'out_synth_working'
    
    # Copy the relevent synth and map files to the files gwas uses
    filename = "synthetic_genomes_" + model_name + ".txt"  
    filepath = data_path / "synthetic_genome_data" / filename
    synth_df = pd.read_csv(filepath, sep=' ')
    newfile = data_path / "synthetic_genomes.txt"
    synth_df.to_csv(newfile, index=False, sep=' ')  
    filename = "map_" + model_name + ".txt"  
    filepath = data_path / "genome_map_data" / filename
    map_df = pd.read_csv(filepath, sep=' ')
    newfile = data_path / "map_batch.txt"
    map_df.to_csv(newfile, index=False, sep=' ')
   
    # Run GWAS
    !rm ../mice_data_set/out_synth/*.csv
    !R --vanilla < ../research_paper_code/src/map_gwas_batch.R &> /tmp/R.log  
    
    filename = data_path / 'pheno_analysis.csv'
    pheno_analysis_df = pd.read_csv(filename)
    pheno_analysis = list(pheno_analysis_df["pheno"])
    
    all_f1s = []

    for phenotype in pheno_analysis:
        
        # Read in the new results

        try:
            synthetic_snps = pd.read_csv(synthetic_gwas_path / f'lm_{phenotype}.csv')  
        except:
            !R --vanilla < ../research_paper_code/src/run_map_batch1.R &> /tmp/R.log
            synthetic_snps = pd.read_csv(synthetic_gwas_path / f'lm_{phenotype}.csv')
        
        synthetic_snps = synthetic_snps.rename(columns={synthetic_snps.columns[0]: 'index'})
        synthetic_snps = synthetic_snps[['index', 'snp', 'p']]
        synthetic_snps['interest'] = synthetic_snps['p'].apply(lambda x: True if x <= 5e-8 else False)
         
        # Read in the original results

        real_snps = pd.read_csv(real_gwas_path / f'lm_{phenotype}_1_79646.csv') #, usecols=['snp', 'p']) # , usecols=['snp', 'p']
        real_snps = real_snps.rename(columns={real_snps.columns[0]: 'index'})
        real_snps = real_snps[['index', 'snp', 'p']]
        real_snps['interest'] = real_snps['p'].apply(lambda x: True if x <= 5e-8 else False)
    

        combined = pd.merge(synthetic_snps, 
             real_snps, 
             how='inner', 
             on=['snp'],
             suffixes=['_synthetic', '_real'])
    
        f1 = round(f1_score(combined['interest_real'], combined['interest_synthetic'], average='weighted'), 4)

        all_f1s.append(f1)
    
    return all_f1s

def save_model_stats(filename):
 
    global models_completed
    global model_info
    global gwas_error_batches
    
    model_names = []
    filenames = []
    sqss = []
    F1s = []
    gwas_errors = []
    
    for model_name in models_complete:
        model_names.append(model_name)
        filenames.append(model_info[model_name]["filename"])
        if "sqs" in model_info[model_name]:
            sqss.append(model_info[model_name]["sqs"])
        else:
            sqss.append(None)
        if "F1" in model_info[model_name]:
            F1s.append(model_info[model_name]["F1"])
        else:
            F1s.append(None)
        if model_name in gwas_error_batches:
            gwas_errors.append(True)
        else:
            gwas_errors.append(False)
            
    results_df = pd.DataFrame({"model_name": model_names, "filename": filenames,
                              "sqs": sqss, "F1": F1s, "gwas_error": gwas_errors})
    
    results_df.to_csv(filename, index=False, header=True)
    print("Completed model stats saved to: " + str(filename))
    
    return results_df


## Start training synthetic models
import time

global pheno_analysis
pass_num = 0

starttime = time.time()

while moreToDo:
    
    pass_num += 1
    print()
    print_pass = "************************************** PASS " + str(pass_num) + " **************************************"
    print(print_pass)
    print("Models training: " + str(models_training_cnt))
    print("Models generating: " + str(models_generating_cnt))
    print("Models complete: " + str(len(models_complete)))
    still_to_start = total_model_cnt - models_training_cnt - models_generating_cnt - len(models_complete)
    print("Models still to start: " + str(still_to_start))
    print("Training errors encountered: " + str(models_error_cnt))
    print("GWAS errors encountered: " + str(gwas_errors))
    print()
    
    # Check for model completion
    chk_model_completion()
    
    # Check for generation completion
    chk_model_generation_completion()
    
    # Check for project completion
    chk_project_completion()
    
    # Start more models if room
    start_more_models()
    
    # Gather complete model stats and save to file
    if len(models_complete) > 0:
        filename_for_stats =  data_path / "Completed_genome_batch_stats.csv"
        results_df = save_model_stats(filename_for_stats)
    
    # Check if we're all done
    if models_complete_cnt == total_model_cnt:
        moreToDo = False
    
    # Sleep for 1  to 5 minutes, adjust to value desired
    time.sleep(60)
    
endtime = time.time() 
exectime = endtime - starttime
exectime_min = round((exectime / 60), 2)
print()
print("************************************** MODELS ALL COMPLETE **************************************")
print(str(total_model_cnt) + " models completed in " + str(exectime_min) + " minutes")
print("A maximum of " + str(MAX_MODELS_RUNNING) + " were allowed to run in parallel")
print(str(models_error_cnt) + " errors occurred which resulted in model retraining")
avg_sec_per_model = exectime / total_model_cnt
avg_min_per_model = round((avg_sec_per_model / 60), 2)
print("Each model took an average of " + str(avg_min_per_model) + " minutes")


# Save all model results
filename_for_stats =  data_path / "synthetic_genome_allbatches_results.csv"
results_df.to_csv(filename_for_stats, index=False, header=True)

# First combine all the synthetic batch results into one dataframe.  This cell can take a while

model_name = models_complete[0]
filename = "synthetic_genomes_" + model_name + ".txt"  
filepath = data_path / "synthetic_genome_data" / filename 
synthetic_genomes = pd.read_csv(filepath, sep=' ')
min_file_length = len(synthetic_genomes.index)

batch_num = 1
while batch_num < len(models_complete):
    model_name = models_complete[batch_num]
    print(batch_num)
    batch_num += 1 
    filename = "synthetic_genomes_" + model_name + ".txt"  
    filepath = data_path / "synthetic_genome_data" / filename 
    synthetic_genomes_batch = pd.read_csv(filepath, sep=' ')
    if (len(synthetic_genomes_batch.index) < min_file_length):
        min_file_length = len(synthetic_genomes_batch.index)
    synthetic_genomes_batch = synthetic_genomes_batch.drop(['id', 'discard'], axis=1)
    synthetic_genomes = pd.concat([synthetic_genomes, synthetic_genomes_batch], axis=1)
            
synthetic_genomes = synthetic_genomes.dropna().reindex
synthetic_genomes= synthetic_genomes.astype({'id': 'int32'})

# Now create of version of map.txt with just the SNPs in this final synthetic dataset
# This will be used in the map code to run GWAS one more time on the overall synthetic dataset
# NOTE: IT'S IMPORTANT TO REMEMBER THE NAME YOU CHOOSE FOR THE MAP FILE.  YOU WILL USE IT IN FURTHER NOTEBOOKS
# I GIVE NO SPECIFIC NAME; IT COULD TAKE THE NAME FROM THE PHENOTYPE TRAIT IN THE ANALYSIS CONFIG FILE.
mapfile = data_path / "map.txt"
mapdata = pd.read_csv(mapfile, sep=' ')

snps = list(synthetic_genomes.columns)
snps.remove("id")
snps.remove("discard")

mapdata_use = mapdata[mapdata["id"].isin(snps)]
filename = "map_allbatches.txt"
mapfile_new = data_path / "genome_map_data" / filename
mapdata_use.to_csv(mapfile_new, sep=' ', header=True, index=False)

# Put the synthetic genome file in chromosome order
map_ids = list(mapdata_use["id"])
col_use = ["id", "discard"]
col_use = col_use + map_ids
synthetic_genomes = synthetic_genomes.filter(col_use)

# Modify the filename to be what you want
# NOTE: IT'S IMPORTANT TO REMEMBER THE NAME YOU CHOOSE FOR THE GENOME FILE.  YOU WILL USE IT IN FURTHER NOTEBOOKS
# I GIVE NO SPECIFIC NAME; IT COULD TAKE THE NAME FROM THE PHENOTYPE TRAIT IN THE ANALYSIS CONFIG FILE.

filename = "synthetic_genomes_allbatches.txt"  
filepath = data_path / "synthetic_genome_data" / filename 
synthetic_genomes = synthetic_genomes.fillna(0)
synthetic_genomes.to_csv(filepath, index=False, header=True, sep=' ')

