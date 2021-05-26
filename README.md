# SMSN

Adapted version of Beaularier's SMSN (See [1]), but with in-sillico control rather than whole-genome amplification, and applied to Sequel I and II data. The SMSN scores are not computed exactly like the ones of Beaulaurier yet - for now it still outputs PacBio's scores.

# Status

Prototype, for research purpose only

- Ran fine on real world data (Sequel I data), Ubuntu 20 LTS x86
-- Not tested in any other platform
- Can be installed and used easily using pip and conda (See "Installation")
- Automated tests are lacking

Since I am a single maintener on this and it was developped for a specific usage on a specific machine, corrections might be required to get it work it on other configurations (e.g Apple computer)

Minor bugs might also occur in specific situations. See "known problems" and "In case of installation problems" for more details

# Software Requirements

- Linux LTS 20 or later (Other might work but are not tested) - x86-64 bits
- conda 4.10

# CLI Usage

 ```console
(smsn) guillaume@A320MA:~$ smsn --help
usage: smsn [-h] --bam BAM --reference REFERENCE
            [--model {SP2-C2,SP3-C3,P6-C4,auto}] --output_csv OUTPUT_CSV
            [--CCS CCS] [--min_identity MIN_IDENTITY]
            [--min_subreads MIN_SUBREADS] [--tmpdir TMPDIR]
            [--verbosity {DEBUG,INFO,WARNING,ERROR,CRITICAL}]
            [--progress_bar PROGRESS_BAR] [--nb_proc NB_PROC]
            [--sizechunks SIZECHUNKS] [--add_context ADD_CONTEXT]
            [--preserve_tmpdir] [--idQvs IDQVS]

optional arguments:
  -h, --help            show this help message and exit
  --bam BAM, -b BAM     Path to a .bam file with all the subreads (adapters
                        sequences must be already removed)
  --reference REFERENCE, -r REFERENCE
                        Path to a genome reference (fasta file).
  --model {SP2-C2,SP3-C3,P6-C4,auto}, -m {SP2-C2,SP3-C3,P6-C4}
                        [REQUIRED] Choose the model for IPD prediction. ]
  --output_csv OUTPUT_CSV, -o OUTPUT_CSV
                        Ouput file (csv) of the methylation analysis. See the
                        README for further details on the output's format.
  --CCS CCS, -c CCS     [FACULTATIVE] Path to the circular consensus
                        corresponding to the .bam subreads. Default = CCS will
                        be recreated from scratch.
  --min_identity MIN_IDENTITY, -i MIN_IDENTITY
                        minimum identity (percentage) of the CCS required to
                        launch analysis on a hole.[DEFAULT : 0.99]. Must be in
                        ]0;1]
  --min_subreads MIN_SUBREADS
                        Minimum number of subreads required to launch analysis
                        on the hole. DEFAULT = 50 (so that its possible to
                        have >=25X per strand on at least one position).
  --tmpdir TMPDIR, -t TMPDIR
                        Tmp directory
  --verbosity {DEBUG,INFO,WARNING,ERROR,CRITICAL}, -v {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        Choose your verbosity. Default: DEBUG
  --progress_bar PROGRESS_BAR, -p PROGRESS_BAR
                        Displays a progress bar. Disabled automatically if
                        verbosity is set to debug[DEFAULT]: False
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
  --add_context ADD_CONTEXT
                        TIn the output .csv file, displays the +12/-12 context
                        around the nucleotide.(Files generated can be
                        sensitively heavier) [DEFAULT: True, possible = True
                        or False]
  --preserve_tmpdir     Forbids deletion of tmp dir (experimental / deprecated
                        / debug only)
  --idQvs IDQVS         Outputs PacBio's identificationQV [DEFAULT: TRUE]
```


**You need to indicate the in-sillico model you'd like to usein parameters. See the dedicated section to known which one you should use**. 


Only SP2-C2 has been tested for now.


#  <a name="whichmodel"></a> Which model should I use ?

As described in the section above, **the --model argument is actually required** in the default configuration. You have to indicate yours

Available models are:

* P6-C4
* SP2-C2
* SP3-C2

**SP2-C2** is recommended for **Sequel I** chemistries
**SP3-C3** is recommanded for **Sequel IIv2** chemistries
**Most RSII user would probably want to use P5-C3 model, which is not supported yet**

> See [here](https://github.com/PacificBiosciences/kineticsTools/pull/71) for more info.


# Installation 

## TL; DR

```console
git clone https://github.com/GDelevoye/SMSN.git
conda env create -n smsn -f ./SMSN/environment.yml
conda activate smsn
pip install -e ./SMSN/
```


## In case of installation problems 

SInce the first release, I was forced to delete my dependencies on conda-forge because it had became prohibitively slow to install smsn.
In case of future problems (installation is too long, solving environment is taking forever, dependencies are broken):

1. The requirements listed in environment.yml must be respected carefully (the more important being **python 3.7**)
2. Don't install smsn in an already existing environment
3. Use [mamba](https://github.com/mamba-org/mamba) if conda is really slow. See  https://github.com/conda/conda/issues/7239 for more info .
4. Delete conda-forge and conda-metachannel from your channels 
5. Use a strict channel priority using  *conda config --set channel_priority false*
6. Update conda through *conda update -n base conda* (Tested with *conda 4.10.1*)
7. In your conda configuration, set the "defaults" conda channel at highest priority, followed by the "bioconda" channel, and remove all the others

Major environment changes were done for 1.0.1
Older commit, however, might help you to find suitable settings that used to work with conda-forge.

# Hardware requirements and calculation time

- SSD is adviced
- 1 recent CPU handles a hole in about ~5 to 10s on average (order of magnitude)
- The size of the files generated varies a lot, but can easily reach several GB per run. 
- min 2GB RAM per processor allocated to the job and 0.13Mb per hole 
-- Each of these two conditions must be met, but don't need to be added



# Known problems

SMSN is a prototype and still has some issues.

* **When SMSN crashes because there's not enough RAM on Linux, don't expect to have an informative log. It will just crash**

* **Multiprocessing deadlocks can occur, especially if --min_subreads < 50 or --min_identity is too low**

    * This is due to PacBio's kineticsTools / ipdSummary not outputing the .csv or .gff file when the coverage is too low

    * This is by far the most annoying problem right now, because the program hangs forever without printing any warning when it happens

    * **To avoid it, I recommand that you don't use --min_subreads < 50 nor --min_identity < 0.99**
 
* **When put in "auto" mode and CCS are not provided by user (default behaviour), SMSN crashes because ipdSummary doesn't recognize the chemistry**

    * Two solutions:

        * Provide consensus built with CCS 3.0.0 manually **or** 
        
        * Specify the model manually with the --model argument (see next section to pick which)

        * The warnings and errors are not implemented correctly
    
* **Bugs can occur when the coverage is really low for all HoleID**

* **Installation with conda can be slow**, especially if you're using conda-forge or conda-metachannel (see the dedicated section to fix it). 

* **Only a very tiny proportion of the code is properly tested for now (prototype)**

* **Some DEBUG/INFO lines might be a bit wrong / imprecise**


# References

[1] Beaulaurier, J., Zhang, X. S., Zhu, S., Sebra, R., Rosenbluh, C., Deikus, G., ... & Fang, G. (2015). Single molecule-level detection and long read-based phasing of epigenetic variations in bacterial methylomes. Nature communications, 6(1), 1-12.
