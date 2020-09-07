# ronaldo
Comparison tool of SARSCOV2 seq and assays

## Install
```
pip install -r requirements.txt 
python ronaldo/ronaldo.py -h 
```

## How to run 
There are two stages to run, first to calculate the required metrics (genome coverage, genome recovery, number of reads \[Illumina\]) and then to filter based on cut-offs. 
If you already have the required statistics you can skip to the filter step. 

### Calculating run metrics 
The calculate module helps generate values for the different metrics. The required inputs are the directory of BAM files for a given run (specified with runname)
and followed by the filename of the Blank in that run (multiple BAM files BLANKS can be specified - the maximum value is retained). 

```
usage: ronaldo.py calculate [-h] [-d DB] [--ctdata CTDATA]
                            [--blank_read_cutoff BLANK_READ_CUTOFF]
                            [--blank_recovery_cutoff BLANK_RECOVERY_CUTOFF]
                            [--ont] [-l READLEN]
                            runname bamfolder N [N ...]

positional arguments:
  runname               Informative label for this run
  bamfolder             Folder of SARSCOV2 BAM files
  N                     Negative control BAM file

optional arguments:
  -h, --help            show this help message and exit
  -d DB, --db DB        DB directory [default: ronaldo_db]
  --ctdata CTDATA       Path to table with assay information
  --blank_read_cutoff BLANK_READ_CUTOFF
                        Run skipped if blanks have number of mapped reads
  --blank_recovery_cutoff BLANK_RECOVERY_CUTOFF
                        Run skipped if blanks have higher perc. genome
                        recovery
  --ont                 Data is OXFORD NANOPORE
  -l READLEN, --readlen READLEN
                        Minimum length for a mapped read (use with ILLUMINA
                        only)
```

The Blanks are used first to check if there is too much hCOV material (e.g. >500 reads or 5% genome recovery by default). If so, no output file will be generated. 
If not, the program will read the remaining BAM files and calculate genome coverage, genome recovery, number of reads \[for Illumina\]). These values will be written
to the specified db directory.

You should also specify a metadata file with the CT information, it should look like this, with the following headers:

| sequencing_platform | sample_name | ct_platform_1  | ct_platform_2  | max_ct_value | min_ct_value | filename              |
|---------------------|-------------|----------------|----------------|--------------|--------------|-----------------------|
| ILLUMINA            | MY_SAMPLE1  | AUSDIAGNOSTICS | AUSDIAGNOSTICS | 15           | 15           | my_sample1.sorted.bam |

Minimum usage should therefore be:

```
python ronaldo/ronaldo.py calculate --ctdata fixed.dat 20200729  Covid-19_Seq/ncovIllumina_sequenceAnalysis_readMapping  BLANK.sorted.bam BLANK2.sorted.bam

```


### Filtering for a final output 

```
FULL USAGE:

usage: ronaldo.py filter [-h] [-d DB] [-o OUTPUT] [-c COVERAGE] [-r RECOVERY]
                         [-n NOREADS] [-t TOTALREADS]
                         sitename

positional arguments:
  sitename              Informative label for your site

optional arguments:
  -h, --help            show this help message and exit
  -d DB, --db DB        DB directory
  -o OUTPUT, --output OUTPUT
                        output directory
  -c COVERAGE, --coverage COVERAGE
                        Minimum fold genome coverage for mapped reads (default: 2)
  -r RECOVERY, --recovery RECOVERY
                        Minimum fold genome recovery for mapped reads (default: 2)
  -n NOREADS, --noreads NOREADS
                        Minimum fold number of mapped reads (default: 5) [Illumina only]
  -t TOTALREADS, --totalreads TOTALREADS
                        Minimum total number of mapped reads (default: 30) [Illumina only]
```

The filtering module is expecting a specified directory (db), with multiple csv files. Each csv file 
could be describing a run (generated from "calculate"). There should be one record/row per sample. There 
must be the following columns with the exact names below:

* **sample_name**	: Unique ID (COG ID)
* **runname**	: A descriptive name of a run. 
* **filename**	: The corresponding BAM filename for this record
* **sequencing_platform**	: ILLUMINA or OXFORD_NANOPORE only 
* **ct_platform_1** : Platform for 1st CT test e.g. AUSDIAGNOSTICS, ROCHE, HOLOLOGIC.
* **ct_platform_2**	: Platform for 2nd CT test e.g. AUSDIAGNOSTICS, ROCHE, HOLOLOGIC.
* **max_ct_value**	: Of the two tests, what is the maximum CT value
* **min_ct_value**	: Of the two tests, what is the minimum CT value
* **blank_coverage**	: mean sequencing coverage for the BLANK for this run. 
* **blank_recovery_10**	: % of bases that were greater than 10X for the BLANK (genome recovery - used for Illumina)
* **blank_recovery_20**	: % of bases that were greater than 20X for the BLANK (genome recovery - used for ONT)
* **blank_reads**	: Number of mapped reads in the BLANK for this run
* **mean_cov**	: mean sequencing coverage
* **pc_pos_gte_20**	: % of bases that were greater than 20X (genome recovery - used for ONT)
* **pc_pos_gte_10**	: % of bases that were greater than 10X (genome recovery - used for Illumina)
* **no_reads** : Number of mapped reads 

The cut-offs are scaled as per the amount detected in the corresponding blank, so it is important to provide 
those values. Adding the run name also helps us detect any batch effects. There are optional flags to change
each cut-off. You also require a site name to be specified, which is used in the output file name. 

The program will add an extra column "false_postive" giving a TRUE/FALSE based on the cut-offs. It will then
write the values to a csv file in the specified output directory. 

Minimum usage should therefore be:

```
python ronaldo/ronaldo.py filter -d ronaldo_db  -o ronaldo_out  NORW 
```
