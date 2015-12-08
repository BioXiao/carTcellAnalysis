# carTcellAnalysis
This repository contains a mini-pipeline for detecting a gene construct (in this case, CAR) in the RNA-seq read data for a sample. The main tool used in the pipeline is Salmon, designed for extremely fast transcriptome mapping and quantification. Salmon is available [here](https://github.com/COMBINE-lab/salmon).

Other important steps in the pipeline mostly involve formatting and configuring reference data, including the sequence of the target gene construct.

### Requirements

The code assumes you have the binary release of Salmon located at `tools/SalmonBeta-0.5.0_OSX-10.10/bin/`; you might need to modify line 8 of `prepData.sh`, line 32 of the `scripts/carGeneSalmon.py`, or line 5 of `testSalmon.py`, if you're running in Unix or with a different version of Salmon.

You'll also need the `pandas` module installed for Python.

### Prepare reference data

```
gffread -w data/sequence/GRCh38_transcripts.fa \
	-g Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	Homo_sapiens.GRCh38.77.gtf
```

[Genome FASTA](ftp://ftp.ensembl.org/pub/release-77/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz) (Homo\_sapiens.GRCh38.dna.primary\_assembly.fa)  
[Gene model GTF](ftp://ftp.ensembl.org/pub/release-77/gtf/homo_sapiens) (Homo_sapiens.GRCh38.77.gtf)  



### Prepare simulated read data for testing
```
$ bash prepSimData.sh data/sequence/carGeneRaw.txt
```

You might get the following error, which I haven't bothered to debug yet. You should be able to just re-run the script once or twice until it works.

```
ValueError: empty range for randrange() (0,0, 0)
```

Generates the following files:

+ `data/sequence/carGeneParts.fasta`
+ `data/simFastqs/simCarParts.fastq`
+ `data/simFastqs/simCarParts_noMulti.fastq`
+ `data/simFastqs/simCarParts_multiOnly.fastq`



### Prepare transcriptome sequence data

```
$ bash prepData.sh data/sequence/carGeneRaw.txt 
```

Generates the following files:

+ `data/sequence/carTranscript.fasta`
+ `data/sequence/GRCh38_transcripts.fasta`



### Run tests

**Note:** these tests are to better understand the behavior of Salmon with different inputs and configurations; they're not actual tests of individual scripts and functions.

```
$ bash testSalmon.sh
```

Generates a number of subfolders and outputs under `data/tests/`. Test results are summarized in the file `data/tests/testSummary.txt`.



### Run pipeline

Once all the necessary data has been prepared, commands for running the main pipeline script (`carGeneSalmon.py`) can be found in `run.sh`. For each set of samples to be processed, the following inputs should be specified:

+ `INPUT_DIR`: folder of FASTQ files (one for each sample)  
+ `OUTPUT_BASE`: folder path and 'tag' for pipeline results  

The `run.sh` script allows you to set up multiple calls of `carGeneSalmon.py` as a batch; to call the script directly, use the following command:

```
$ python scripts/carGeneSalmon.py $INPUT_DIR $OUTPUT_BASE
```

This script does a two-pass search for the CAR gene for each sample. First, reads are mapped only to the CAR transcript sequence. As mapping is done based on k-mers, this step could lead to false positives from endogenous transcripts, and is thus considered 'low confidence' evidence of CAR gene presense. Second, reads are mapped to an index for the full human (GRCh38) transcriptome, with the CAR transcript added on. In this case, if reads do not map *uniquely* to CAR (e.g., segment-spanning reads or reads mapping to non-human segments), the Salmon model should estimate a lower or zero abundance for the CAR transcript.

Generates the following files:

+ `data/results/<tag>.csv`
+ `data/results/<tag>_high_conf.csv`


### Inspect results

Some example R code for merging and summarizing results can be found in `inspectResults.R`.