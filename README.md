# Interface for Gretel AI

The purpose of this repo is to create an interface to use [GretelAI](https://gretel.ai/) for synthetic data generation. 
The pipeline collects the source code from [gretelai/synthetic data genomics](https://github.com/gretelai/synthetic-data-genomics) to create and build phenomes and genomes, uses mouse example data, and allows the user to extend to choose phenotype traits and covariates of interest, as well as model behaviour.

The following shows how to use it on a test case.

## Installation
Clone the repository; ceate and activate a conda environment:
```
conda env create -n gretel -f env/gretel_env.yaml
conda activate gretel
```


## Download example data
Download mouse example data:
```
./config/download.sh
```

## Prepare configuration files
Two files are needed: 'config/phenome_analysis.yml' specifies the phenotypes and covariates you want to analyse; and 'config/model_params.yml' specifies the model parameters. The names and location of these files should not be modified.


The phenome_analysis.yml config file in this test case focuses on traits and covariates studied in this [paper](https://doi.org/10.1038/ng.3609). The user could choose between several traits: 'TA', 'SW16', 'tibia', 'EDL', 'soleus', 'plantaris', 'gastroc', 'SW6', 'sacweight', 'BMD', 'abBMD'.

## Run pipeline
Run whole pipeline with:
```
./src/gretel.sh
```

or run each module individually:
```
python src/01_create_synthetic_phenomes.py
python src/02_build_synthetic_phenomes.py
python src/03_create_synthetic_genomes.py
python src/04_build_synthetic_genomes.py
```

# Remarks
The goal of this repo was to create an easy-to-use interface for gretel-ai. My idea was to create such interface as a sequential call of the different steps for creating synthetic phenome and genome data. The user would be able to supply phenotype and model configuration parameters through json files. The pipeline itself is already written in several study cases provided by gretel-ai. Unfortunately I could not test the idea, as I run out of the free credits needed to run gretel. 