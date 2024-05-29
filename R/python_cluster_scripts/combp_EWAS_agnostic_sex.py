# # In bash, set path where comb-p is located
# conda activate DMR
# export PYTHONPATH=~/.local/lib/python3.8/site-packages/cpv-0.50.3-py3.8.egg/cpv/
# # Check PYTONPATH
# echo $PYTHONPATH

# # In bash, run python sepages/phenols_phthalates_DNAm/R/combp_EWAS_agnostic_sex.py (code below)

# Import modules
import os
import os.path

# Set path where .BED files are located
path = "/home/slavakp/sepages/phenols_phthalates_DNAm/analysis/agnostic_sex/bed_files"
files = os.listdir(path = path)
# For each .BED file in the path, run comb-p pipeline analysis
for bed_file in files:
    path_to_bed = os.path.join(path, bed_file)
    prefix = bed_file.split('.')[0] + ""
    dmr = "comb-p pipeline -c 4 --seed 0.001 --dist 500 -p " + prefix + " --region-filter-p 0.05 --region-filter-n 2 " + path_to_bed
    os.system(dmr)
