# THIS WRAPPER CALLS EACH MODULE SEQUENTIALLY.

# Get main folder
DIR="$(cd "$(dirname "$0")" && pwd)"

# PATH TO EACH MODULE OF THE PIPELINE.
CREATE_PHENOMES=${DIR}/src/01_create_synthetic_phenomes.py
BUILD_PHENOMES=${DIR}/src/02_build_synthetic_phenomes.py
CREATE_GENOMES=${DIR}/src/03_create_synthetic_genomes.py
BUILD_GENOMES=${DIR}/src/04_build_synthetic_genomes.py
# COMPARE=${DIR}/src/05_compare_associations.py

# ------------------------------------------------------------------- #

echo "Part 1 - Creating synthetic phenomes..."
python ${CREATE_PHENOMES}
echo "#-------------------------------------#"
echo "Part 2 - Building synthetic phenomes..."
python ${BUILD_PHENOMES}
echo "#-------------------------------------#"
echo "Next steps were untested, as free gretel credits were used up. Stop here."
# echo "Part 3 - Building synthetic genomes..."
# python ${CREATE_GENOMES}
# echo "#-------------------------------------#"
# echo "Part 4 - Building synthetic genomes..."
# python ${BUILD_GENOMES}
