# SMSN

Version of Beaularier's SMSN, but with in-sillico control rather than whole-genome amplification, and applied to Sequel I and II data

# Status

Pre-alpha - Seems to work but not tested

# Software Requirements

- Linux LTS 20 or later (Other might work but are not tested) - x86-64 bits
- conda

# Hardware requirements and computational times

- Count about 15.2*2 = 30.4 Mb par hole must be free on the hard drive
- SSD is adviced, since lots of things happen on the hard drive
- On a ryzen 5 2600 + HDD 7200tr/mn, 1 CPU handles a hole in about ~10s on average
- min 2GB RAM per processor allocated to the job and 0.13Mb per hole (Each of these two conditions must be met, but don't need to be added)

# Installation 

Because the program relies on ***very precise*** versions of PacBio's tools, python 3.7 **MUST** be used and all the requirements listed in environment.yml ***must*** be respected carefully. The only easy way of not getting wrong is by using a virtual environment manager - e.g conda.

```console
git clone https://github.com/GDelevoye/SMSN.git
conda env create -n smsn -f ./SMSN/environment.yml
conda activate smsn
pip install -e ./SMSN/
```

Known issue: Conda takes LOTS of time to build everything. This is due to the conda solver, and shoudl be solved with the release of conda 5.0. However, even if it's slow, it works well after ~ 1 hour of installation on test machines. See  https://github.com/ContinuumIO/anaconda-issues/issues/9480 for more info . Building the environment with https://github.com/mamba-org/mamba might help if installing is really too slow.

# Usage

```console
(smsn) guillaume@A320MA:~$ smsn --help
usage: smsn [-h] --bam BAM [--CCS CCS] --reference REFERENCE --output_csv
            OUTPUT_CSV [--min_identity MIN_IDENTITY]
            [--min_subreads MIN_SUBREADS] [--tmpdir TMPDIR]
            [--verbosity {DEBUG,INFO,WARNING,ERROR,CRITICAL}] [--progress_bar]
            [--nb_proc NB_PROC] [--sizechunks SIZECHUNKS]

optional arguments:
  -h, --help            show this help message and exit
  --bam BAM, -b BAM     Path to a .bam file with all the subreads (adapters
                        sequences must be already removed)
  --CCS CCS, -c CCS     [FACULTATIVE] Path to the circular consensus
                        corresponding to the .bam subreads. Default = CCS will
                        be recreated from scratch.
  --reference REFERENCE, -r REFERENCE
                        Path to a genome reference (fasta file).
  --output_csv OUTPUT_CSV, -o OUTPUT_CSV
                        Ouput file (csv) of the methylation analysis. See the
                        README for further details on the output's format.
  --min_identity MIN_IDENTITY, -i MIN_IDENTITY
                        minimum identity (percentage) of the CCS required to
                        launch analysis on a hole.
  --min_subreads MIN_SUBREADS, -s MIN_SUBREADS
                        Minimum number of subreads required to launch analysis
                        on the hole. DEFAULT = 50 (so that its possible to
                        have >=25X per strand on at least one position)
  --tmpdir TMPDIR, -t TMPDIR
                        Tmp directory
  --verbosity {DEBUG,INFO,WARNING,ERROR,CRITICAL}, -v {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        Choose your verbosity. Default: DEBUG
  --progress_bar, -p    Displays a progress bar. Disabled automatically if
                        verbosity is set to debug[DEFAULT]: Set to True if
                        DEBUG as verbosity, otherwise FALSE
  --nb_proc NB_PROC, -n NB_PROC
                        Multiprocessing on n CPU. Default: 1.
  --sizechunks SIZECHUNKS, -k SIZECHUNKS
                        Because the subreads .bam will often not fit entirely
                        in RAM and because the methylation analysis itself
                        generates lots of data, SMSN will pause and perform
                        intensive I/O operation every S holes it has analyzed.
                        Lower values are better for machines that are limited
                        in RAM. The optimal nb_proc/sizehunks ratio will vary
                        from one computer to another. In case sizechunks <
                        nb_proc, SMSN will use sizechunks = 20x nb_proc
                        instead.
```
