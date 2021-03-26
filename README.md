# SMSN

Version of Beaularier's SMSN, but with in-sillico control rather than whole-genome amplification, and applied to Sequel I and II data

# Installation 

Because the program relies on precise versions of PacBio's tools, python 3.7 MUST be used - preferably within a conda environment.

```console
git clone https://github.com/GDelevoye/SMSN.git
conda env create -n smsn -f ./SMSN/environment.yml
conda activate smsn
pip install -e ./SMSN/
```

# Known issue during installation

You might sometimes encounter this:

```console
(base) [bioclusts01 ~ 10:57]$ conda create -n smsn -f GitHub/SMSN/environment.yml 
Solving environment: failed

PackagesNotFoundError: The following packages are not available from current channels:

  - github/smsn/environment.yml

Current channels:

  - https://conda.anaconda.org/bioconda/linux-64
  - https://conda.anaconda.org/bioconda/noarch
  - https://conda.anaconda.org/r/linux-64
  - https://conda.anaconda.org/r/noarch
  - https://repo.anaconda.com/pkgs/main/linux-64
  - https://repo.anaconda.com/pkgs/main/noarch
  - https://repo.anaconda.com/pkgs/free/linux-64
  - https://repo.anaconda.com/pkgs/free/noarch
  - https://repo.anaconda.com/pkgs/r/linux-64
  - https://repo.anaconda.com/pkgs/r/noarch
  - https://repo.anaconda.com/pkgs/pro/linux-64
  - https://repo.anaconda.com/pkgs/pro/noarch
  - https://conda.anaconda.org/conda-forge/linux-64
  - https://conda.anaconda.org/conda-forge/noarch
```

Whis is a conda bug, that makes it ignore the channels in some conda versions. You can either:
- Update conda to a newer version
- Add manually the required channels on your .condarc file :

```console
conda config --add channels default conda-forge bioconda
```


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
