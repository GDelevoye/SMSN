# SMSN

Adapted version of Beaularier's SMSN (See [1]), but with in-sillico control rather than whole-genome amplification, and applied to Sequel I and II data. The SMSN scores are not computed exactly like the ones of Beaulaurier yet - for now it still outputs PacBio's scores.

# Status

Pre-alpha / Prototype

- Seems to work correctly with SequelI data
- Lots of tests are lacking

# /!\ IMPORTANT /!\ Known problems

SMSN is a prototype and still has some issues.

**- Multiprocessing deadlocks can occur, especially if --min_subreads < 50 or --min_identity is too low**

-- This is due to PacBio's kineticsTools / ipdSummary not outputing the .csv or .gff file when the coverage is too low
-- This is by far the most annoying problem right now, because the program hangs forever without printing any warning when it happens
-- **To avoid it, I recommand that you don't use --min_subreads < 50 nor --min_identity < 0.99**

- **When put in "auto" mode and CCS are not provided by user (default behaviour), SMSN crashes because ipdSummary doesn't recognize the chemistry**

-- Either provide consensus built with CCS 3.0.0 manually **or** specify the model manually (see next section)
-- The warnings and errors are not implemented correctly

- Some DEBUG/INFO lines might be a bit wrong / imprecise yet

- Only a very tiny proportion of the code is properly tested for now (prototype)

- ipdSummary doesn't output anything for low coverage. If no holeID has sufficiently high effective coverage in the input .bam file, the program might close properly with a decent explicative warning/error. But it might also crash savagely (I didn't tested it sufficiently to know)

- Installation with conda is super-slow (see the dedicated section). It relies on some specific commit of GitHub repositories. 
 
--If they are unreachable (i.e: GitHub is down, the repository is set in private), which should never happen (but we never know) then the installation would fail and dependencies would have to be installed manually. Everything runs fine on the 03/31/2021. 

#  <a name="whichmodel"></a> Which model should I use ?

Available models are:

* P6-C4
* SP2-C2
* SP3-C2

**SP2-C2** is recommended for **Sequel I** chemistries
**SP3-C3** is recommanded for **Sequel IIv2** chemistries
**Most RSII user would probably want to use P5-C3 model, which is not supported yet**

> See [here](https://github.com/PacificBiosciences/kineticsTools/pull/71) for more info.


# Software Requirements

- Linux LTS 20 or later (Other might work but are not tested) - x86-64 bits
- conda 4.7.10 or later

# Hardware requirements and calculation time

- Count about 15.2Mb par hole must be free on the hard drive in both the tmp_dir and the directory where the .csv output must be produced
- SSD is adviced, since lots of things happen on the hard drive
- On a ryzen 5 2600 + HDD 7200tr/mn, 1 CPU handles a hole in about ~5 to 10s on average depending on the parameters and input files
- The size of the files generated varies a lot, but can easily reach several GB per run. Make sure you have the free space
- min 2GB RAM per processor allocated to the job and 0.13Mb per hole 
-- Each of these two conditions must be met, but don't need to be added

**WARNING: When SMSN crashes because there's not enough RAM on Linux, don't expect to have an informative log. It will just crash**

# Installation 

Because the program relies on ***very precise*** versions of PacBio's tools, python 3.7 **MUST** be used and all the requirements listed in environment.yml ***must*** be respected carefully. The only easy way of not getting wrong is by using a virtual environment manager - e.g conda.

```console
git clone https://github.com/GDelevoye/SMSN.git
conda env create -n smsn -f ./SMSN/environment.yml
conda activate smsn
pip install -e ./SMSN/
```

Known issue: Conda takes LOTS of time to build everything. This is due to the conda solver, and shoudl be solved with the release of conda 5.0. However, even if it's slow, it works well after ~ 1 hour of installation on test machines. See  https://github.com/conda/conda/issues/7239 for more info . Building the environment with https://github.com/mamba-org/mamba might help if installing is really too slow.

# Usage

Only two things are required: 

* A .bam with your PacBio Sequel I or II subreads (where adapter sequences have been removed)
* A .fasta reference of your genome of interest. 


The PacBio tools parse automatically the headers to use the right in-sillico models.

 ```console
(smsn) guillaume@A320MA:~$ smsn --help
usage: smsn [-h] --bam BAM [--CCS CCS] --reference REFERENCE --output_csv
            OUTPUT_CSV [--min_identity MIN_IDENTITY]
            [--min_subreads MIN_SUBREADS] [--tmpdir TMPDIR]
            [--verbosity {DEBUG,INFO,WARNING,ERROR,CRITICAL}] [--progress_bar]
            [--nb_proc NB_PROC] [--sizechunks SIZECHUNKS] [--preserve_tmpdir]

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
  --preserve_tmpdir     Forbids deletion of tmp dir (experimental / deprecated
                        / debug only)
```

# References

[1] Beaulaurier, J., Zhang, X. S., Zhu, S., Sebra, R., Rosenbluh, C., Deikus, G., ... & Fang, G. (2015). Single molecule-level detection and long read-based phasing of epigenetic variations in bacterial methylomes. Nature communications, 6(1), 1-12.
