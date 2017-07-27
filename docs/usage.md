Usage
-----

### Display Help Message

The `help` parameter displays the available options and commands.
```
$ nextflow run auir.nf --help
```

### File Inputs

#### Set custom sequence data

The `reads` parameter accepts sequence files in standard fastq and gz format.
```
$ nextflow run auir.nf --reads "data/raw/*_R{1,2}.fastq"
```

#### Set host genome

The `host` parameter accepts a fasta formatted host genome.
```
$ nextflow run auir.nf --reads "data/raw/*_R{1,2}.fastq" --host "data/host/gallus.fa"
```

#### Set host index

The `index` parameter allows you to upload pre-built host indexes produced by BWA.
```
$ nextflow run auir.nf --reads "data/raw/*_R{1,2}.fastq" --host "data/host/gallus.fa" --index "data/index/*"
```

### File Outputs

#### Set output directory

The `output` parameter writes output files to the desired directory.
```
$ nextflow run auir.nf --reads "data/raw/*_R{1,2}.fastq" --output "test"
```

### Other Parameters

#### Set custom trimming parameters

```
$ nextflow run auir.nf \
    --reads "data/raw/*_R{1,2}.fastq" \
    --leading 3 \
    --trailing 3 \
    --minlen 36 \
    --slidingwindow 4 \
    --adapters "data/adapters/nextera.fa" \
    --output "test"
```

#### Set number of threads

The `threads` parameter specifies the number of threads to use for each process.
```
$ nextflow run auir.nf --reads "data/raw/*_R{1,2}.fastq" --threads 16 --output "test"
```
