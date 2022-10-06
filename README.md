[![CircleCI](https://circleci.com/gh/EMeyerLab/SMSN/tree/clean_article.svg?style=shield&circle-token=bfc806f868edcd05f6f15e1041c290ea11f24531)](https://app.circleci.com/pipelines/github/EMeyerLab/SMSN)

# SMSN

Adapted version of Beaularier's SMSN

See the reference [1] and [this Github repo](https://github.com/fanglab/SMALR). The differences are the following :

* Here, we use the PacBio in-sillico control rather than a whole-genome amplified DNA.
* We applied our pipeline to Sequel I, rather than RS II
* Not tested on Sequel II data but should work too.
* The SMSN scores are not computed exactly like the ones of Beaulaurier et al.

Our [companion article](https://github.com/GDelevoye/SMSN/blob/main/article/preprint/SMSN_Article_preprint_v01.pdf) explains how the scores are computed and how they can be interpreted

# Pre-print

A preprint describing the software's architecture, goals and results in _E. coli_ is available [here](https://github.com/GDelevoye/SMSN/blob/main/article/preprint/SMSN_Article_v2.pdf). Please read [the disclaimer](https://github.com/GDelevoye/SMSN/blob/main/article/preprint/DISCLAIMER.md) as well.


# Software Requirements

- Linux LTS 20 or later (Other might work but are not tested) - x86-64 bits
- conda >= 4.10.3

# Installation (TL;DR)


```console
git clone https://github.com/GDelevoye/SMSN.git
conda activate
conda create -n smsn
conda env update -n smsn -f ./SMSN/environment.yml
conda activate smsn
pip install -e ./SMSN/
```

## In case of installation problems 

 * Execute the installation script line by line in a terminal rather than running it in a text file.
 * Installation can sometimes be slow due to conda-forge. See  [the detailed explanations](#in-case-of-installation-problems) for more info. Workardounds :
 * Don't install smsn in an already existing environment !
 * Use [mamba](https://github.com/mamba-org/mamba) if conda is really slow. See  https://github.com/conda/conda/issues/7239 for more info .
 * Delete conda-forge and conda-metachannel from your channels if conda is too slow
 * Use a strict channel priority using  *conda config --set channel_priority false* if conda is too slow
 * Update conda through *conda update -n base conda* (Tested with *conda 4.10.1*)
 * In your conda configuration, set the "defaults" conda channel at highest priority, followed by the "bioconda" channel, and remove all the others if conda is too slow

# CLI Usage


**Mandatory inputs**:

1. A .bam file with the pacbio subreads
2. A .fasta file of your reference genome
3. A model name for the in-sillico control (See next section)
4. The path for an output .csv file


 ```console
(smsn) user@computer:~$ smsn --help
usage: smsn [-h] --bam BAM --reference REFERENCE --model {SP2-C2,SP3-C3,P6-C4}
            --output_csv OUTPUT_CSV [--CCS CCS] [--min_identity MIN_IDENTITY]
            [--min_subreads MIN_SUBREADS] [--tmpdir TMPDIR]
            [--verbosity {DEBUG,INFO,WARNING,ERROR,CRITICAL}]
            [--progress_bar PROGRESS_BAR] [--nb_proc NB_PROC]
            [--sizechunks SIZECHUNKS] [--add_context ADD_CONTEXT]
            [--idQvs IDQVS]

This software implements the SMSN approach of PacBio sequencing in the same
way as Beaulaurier et al 2015, but contrary to Beaulaurier's software, this
one can be used even in the absence of a PCR-amplified control, and on more
recent data (tested on Sequel I).

DELEVOYE Guillaume 2020.
thttps://github.com/EMeyerLab/SMSN

required arguments:
  --bam BAM, -b BAM     Path to a .bam file with all the subreads (adapters
                        sequences must be already removed)
  --reference REFERENCE, -r REFERENCE
                        Path to a genome reference (fasta file).
  --model {SP2-C2,SP3-C3,P6-C4}, -m {SP2-C2,SP3-C3,P6-C4}
                        Choose the model for IPD prediction. [DEFAULT: auto
                        (PacBio's kineticsTools softwarechoses after it has
                        parsed the input file)]
  --output_csv OUTPUT_CSV, -o OUTPUT_CSV
                        Ouput file (csv) of the methylation analysis. See the
                        README for further details on the output's format.

optional arguments:
  --CCS CCS, -c CCS     [FACULTATIVE] Path to the circular consensus
                        corresponding to the .bam subreads. Default = CCS will
                        be recreated from the subreads provided.
  --min_identity MIN_IDENTITY, -i MIN_IDENTITY
                        minimum identity (percentage) of the CCS required to
                        launch analysis on a hole.[DEFAULT : 0.99]. Must be in
                        ]0;1]
  --min_subreads MIN_SUBREADS
                        Minimum number of subreads required to launch analysis
                        on the hole. DEFAULT = 50 (so that its possible to
                        have >=25X per strand on at least one position).
  --tmpdir TMPDIR, -t TMPDIR
                        Tmp directory (DEFAULT : smsn_tmpdir_[DATE_HOUR])
  --verbosity {DEBUG,INFO,WARNING,ERROR,CRITICAL}, -v {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        Choose your verbosity. Default: INFO
  --progress_bar PROGRESS_BAR, -p PROGRESS_BAR
                        Displays a progress bar. Disabled automatically if
                        verbosity is set to debug[DEFAULT]: False
  --nb_proc NB_PROC, -n NB_PROC
                        Multiprocessing on n CPU. Default: 1.
  --sizechunks SIZECHUNKS, -k SIZECHUNKS
                        The subreads will often not fit entirely in RAM and
                        the methylation analysis itself generates lots of I/O
                        usage. Because of it, smsn will pause and perform
                        intensive I/O operation every S holes it has analyzed.
                        Lower values are better for machines that are limited
                        in RAM. The optimal nb_proc/sizehunks ratio will vary
                        from one computer to another. In case sizechunks <
                        nb_proc, SMSN will use sizechunks = 20x nb_proc
                        instead. DEFAULT : 5000
  --add_context ADD_CONTEXT
                        In the output .csv file, displays the +12/-12 context
                        around the nucleotide.(Files generated can be
                        sensitively heavier) [DEFAULT: True, choices = True or
                        False]
  --idQvs IDQVS         Outputs PacBio's identificationQV [DEFAULT: True,
                        choices = True or False]

```

# Development

## Roadmap 🛠️

Prototype, for research purpose only (see License). 

Here are the possible additional developments :

- Better documentation about the different in-silico controls
- conda package release 
- Check the pipeline with Older (RS II) or newer (Sequel II.v2) data

## Tests 

There are two ways to run the test : 1) Either you want to test the repository or 2) You want to test your pip installation's correctness. 

### Test your local copy of the repository

```console 
pytest . # Docstrings should be tested too because of the pytest.ini file at the root
```

To run the same tests with a coverage analysis :

```console
coverage run -m pytest . -vv && coverage report && coverage html
```

### Test the pip installation

Testing files are not shipped with the pip installation. You can however test your pip-installation from the repository. If you have made a non-editable pip installation, testing your pip installation can be done running : 

```console
coverage run -m pytest tests/ --doctest-modules --pyargs smsn -vvv && coverage report && coverage html
```

# FAQ 

I was the only user of this software, so this is not really a FAQ, but rather an anticipation of what problems you might face when using this code.

##  <a name="whichmodel"></a> Which model should I use ?

As described in the section above, **the --model argument is actually required** in the default configuration. You have to indicate yours

Available models are:

* P6-C4
* SP2-C2
* SP3-C2

**SP2-C2** is recommended for **Sequel I** chemistries
**SP3-C2** is recommanded for **Sequel IIv2** chemistries
**Most RSII user would probably want to use P5-C3 model, which is not supported**

Only SP2-C2 (Sequel I) has been tested yet, but the others should work too.

> See [here](https://github.com/PacificBiosciences/kineticsTools/pull/71) for more info.

# References

[1] Beaulaurier, J., Zhang, X. S., Zhu, S., Sebra, R., Rosenbluh, C., Deikus, G., ... & Fang, G. (2015). Single molecule-level detection and long read-based phasing of epigenetic variations in bacterial methylomes. Nature communications, 6(1), 1-12.
