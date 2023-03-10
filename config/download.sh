# Following links are from https://datadryad.org/stash/dataset/doi:10.5061/dryad.2rs41
WORKING_DIR=mice_data_set

# Download experiment datasets
mkdir -p $WORKING_DIR/data $WORKING_DIR/out  $WORKING_DIR/data $WORKING_DIR/out_synth $WORKING_DIR/out_synth_working $WORKING_DIR/out_synth/manh_plots $WORKING_DIR/gemma/output
mkdir -p $WORKING_DIR/data/genome_training_data
mkdir -p $WORKING_DIR/data/genome_map_data
mkdir -p $WORKING_DIR/data/synthetic_genome_data
wget -O $WORKING_DIR/data/geno.txt.gz https://datadryad.org/stash/downloads/file_stream/4344
gunzip -f $WORKING_DIR/data/geno.txt.gz
wget -O $WORKING_DIR/data/map.txt https://datadryad.org/stash/downloads/file_stream/4342
wget -O $WORKING_DIR/data/pheno.csv https://datadryad.org/stash/downloads/file_stream/4340

# Download experiment results
wget -c https://gretel-public-website.s3.amazonaws.com/synthetics/genomics-experiment-output.tar.gz -O - | tar -xz
rm -f genomics-experiment-output.tar.gz

# Download GEMMA
wget https://github.com/genetics-statistics/GEMMA/releases/download/v0.98.4/gemma-0.98.4-linux-static-AMD64.gz -P $WORKING_DIR/gemma/
gunzip -f $WORKING_DIR/gemma/gemma-0.98.4-linux-static-AMD64.gz
chmod u+x $WORKING_DIR/gemma/gemma-0.98.4-linux-static-AMD64

