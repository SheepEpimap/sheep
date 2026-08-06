# 1) Create the environment and install dependencies (mamba is faster; conda works as well)
mamba create -n gat -c conda-forge -c bioconda python=2.7 gat=1.3.6
# 2) Activate the environment
conda activate gat
# 3) Verify the installation
gat-run.py --help
