# carTcellAnalysis

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