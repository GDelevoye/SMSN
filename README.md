# SMSN
Version of Beaularier SMSN with in-sillico control

# Installation 

En python 3.7 on peut trouver tout le nécessaire sur anaconda/bioconda
Dans les autres versions de python, non

```console
conda create -n smsn python=3.7
conda activate smsn
conda install -c bioconda pbcore
conda install -c bioconda pbcoretools
conda install -c bioconda pbcommand
conda install -c bioconda pbbam
conda install -c bioconda blasr
conda install -c bioconda samtools=1.9 --force-reinstall
conda install -c bioconda pbccs
conda install -c bioconda pybigwig
conda install psutil
conda install pandas
pip install pandarallel
git clone https://github.com/GDelevoye/kineticsTools.git
pip install ./kineticsTools/
git clone https://github.com/GDelevoye/SMSN.git
pip install ./SMSN/
```

# Usage

```console
usage: smsn [-h] --bam BAM [--CCS CCS] --reference REFERENCE [--model {SP2-C2,C2,P4-C2,P5-C3,P6-C4,XL-C2,XL-XL,auto}] --output_csv OUTPUT_CSV [--min_identity MIN_IDENTITY] [--min_subreads MIN_SUBREADS]
            [--tmpdir TMPDIR] [--verbosity {DEBUG,INFO,WARNING,ERROR,CRITICAL}] [--progress_bar] [--nb_proc NB_PROC] [--sizechunks SIZECHUNKS]

optional arguments:
  -h, --help            show this help message and exit
  --bam BAM, -b BAM     Path to a .bam file with all the subreads (adapters sequences must be already removed)
  --CCS CCS, -c CCS     [FACULTATIVE] Path to the circular consensus corresponding to the .bam subreads. Default = CCS will be recreated from scratch.
  --reference REFERENCE, -r REFERENCE
                        Path to a genome reference (fasta file).
  --model {SP2-C2,C2,P4-C2,P5-C3,P6-C4,XL-C2,XL-XL,auto}, -m {SP2-C2,C2,P4-C2,P5-C3,P6-C4,XL-C2,XL-XL,auto}
                        Choose the model for IPD prediction. See https://github.com/GDelevoye/ipdtools#-which-model-should-i-use- for more info. DEFAULT: Automatic guess according to the input .bam
  --output_csv OUTPUT_CSV, -o OUTPUT_CSV
                        Ouput file (csv) of the methylation analysis. See the README for further details on the output's format.
  --min_identity MIN_IDENTITY, -i MIN_IDENTITY
                        minimum identity (percentage) of the CCS required to launch analysis on a hole.
  --min_subreads MIN_SUBREADS, -s MIN_SUBREADS
                        Minimum number of subreads required to launch analysis on the hole. DEFAULT = 50 (so that its possible to have >=25X per strand on at least one position)
  --tmpdir TMPDIR, -t TMPDIR
                        Tmp directory
  --verbosity {DEBUG,INFO,WARNING,ERROR,CRITICAL}, -v {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        Choose your verbosity. Default: DEBUG
  --progress_bar, -p    Displays a progress bar. Disabled automatically if verbosity is set to debug
  --nb_proc NB_PROC, -n NB_PROC
                        Multiprocessing on n CPU. Default: 1.
  --sizechunks SIZECHUNKS, -k SIZECHUNKS
                        Because the subreads .bam will often not fit entirely in RAM and because the methylation analysis itself generates lots of data, SMSN will pause and perform intensive I/O operation
                        every S holes it has analyzed. Lower values are better for machines that are limited in RAM. The optimal nb_proc/sizehunks ratio will vary from one computer to another. In case
                        sizechunks < nb_proc, SMSN will use sizechunks = 20x nb_proc instead.

```
