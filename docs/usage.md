Usage
-----

### Display Help Message

```
$ nextflow run auir.nf --help
```

### File Inputs

#### Use custom sequence data
```
$ nextflow run auir.nf --reads "data/raw/*_R{1,2}.fastq"
```

#### Use host genome
```
$ nextflow run auir.nf --reads "data/raw/*_R{1,2}.fastq" --host "data/host/gallus.fa"
```

#### Use host index
```
$ nextflow run auir.nf --reads "data/raw/*_R{1,2}.fastq" --host "data/host/gallus.fa" --index "data/index/*"
```

### File Outputs
```
$ nextflow run auir.nf --reads "data/raw/*_R{1,2}.fastq" --output "test"
```

### Other Parameters
```
$ nextflow run auir.nf \
    --reads "data/raw/*_R{1,2}.fastq" \
    --leading 3 \
    --trailing 3 \
    --minlen 36 \
    --slidiningwindow 4 \
    --adapters "data/adapters/nextera.fa"
    --output "test"
```
