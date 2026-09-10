# System requirements

## Demonstration workflow

The repository demo was tested with Python 3.12.3 on Windows 11. Continuous
integration also runs it on the current GitHub-hosted Ubuntu image with Python
3.11. It requires one CPU, less than 250 MB of RAM, less than 10 MB of temporary
disk space and normally completes in under two seconds after dependencies are
loaded. No GPU, scheduler or network connection is required at run time.

## Full analysis workflows

The complete sequencing and epigenomic analyses were developed for a Linux HPC
environment. Individual scripts use multiple CPU cores and, for large alignment
or genome-wide analyses, tens to hundreds of gigabytes of RAM and substantial
temporary storage. ChromBPNet model training requires a CUDA-capable GPU.
Array-oriented scripts use Slurm (`sbatch`) or site-specific cluster commands;
adapt scheduler directives and resource requests for the target system.

Storage requirements depend on the selected workflow and input datasets. The
complete 516-dataset analysis is not intended as a laptop-scale example. Use
the small demo to verify installation, then review the README and resource
settings in the specific pipeline or figure directory before running real data.
