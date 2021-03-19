# SMSN
Version of Beaularier SMSN with in-sillico control

# Installation 

```console
user@computer:$ git clone https://github.com/GDelevoye/SMSN.git
user@computer:$ pip install ./SMSN
```

# Usage

```console
guillaume@A320MA:~/GitHub$ smsn --help
usage: smsn [-h] --bam BAM --reference REFERENCE [--model {SP2-C2,C2,P4-C2,P5-C3,P6-C4,XL-C2,XL-XL,auto}] [--frequency FREQUENCY] --output_csv OUTPUT_CSV [--tmpdir TMPDIR]
            [--verbosity {DEBUG,INFO,WARNING,ERROR,CRITICAL}] [--progress_bar] [--n_proc N_PROC]

optional arguments:
  -h, --help            show this help message and exit
  --bam BAM, -b BAM     Path to a .bam file with all the subreads (adapters sequences must be already removed)
  --reference REFERENCE, -r REFERENCE
                        Path to a genome reference (fasta file).
  --model {SP2-C2,C2,P4-C2,P5-C3,P6-C4,XL-C2,XL-XL,auto}, -m {SP2-C2,C2,P4-C2,P5-C3,P6-C4,XL-C2,XL-XL,auto}
                        Choose the model for IPD prediction. See https://github.com/GDelevoye/ipdtools#-which-model-should-i-use- for more info. DEFAULT: Automatic guess according to the input .bam
  --frequency FREQUENCY, -f FREQUENCY
                        [NOT MANDATORY] Frequency of the sequencer (Hz). Default = AUTO. In the tested datasets, Sequel I had 80Hz and RSII had 75Hz.
  --output_csv OUTPUT_CSV, -o OUTPUT_CSV
                        Ouput file (csv) of the methylation analysis. See the README for further details on the output's format.
  --tmpdir TMPDIR, -t TMPDIR
                        Tmp directory
  --verbosity {DEBUG,INFO,WARNING,ERROR,CRITICAL}, -v {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        Choose your verbosity. Default: INFO
  --progress_bar, -p    Displays a progress bar
  --n_proc N_PROC, -n N_PROC
                        Multiprocessing on n CPU. Default: 1. If set to -1, will use all available cores.
```
