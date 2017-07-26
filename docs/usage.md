Usage
-----

### Display Help Message

```
$ nextflow run auir.nf --help
```

### File Inputs

#### Set custom sequence data
```
$ nextflow run auir.nf --reads "data/raw/*_R{1,2}.fastq"
```

#### Set host genome
```
$ nextflow run auir.nf --reads "data/raw/*_R{1,2}.fastq" --host "data/host/gallus.fa"
```

#### Set host index
```
$ nextflow run auir.nf --reads "data/raw/*_R{1,2}.fastq" --host "data/host/gallus.fa" --index "data/index/*"
```

### File Outputs

#### Set output directory
```
$ nextflow run auir.nf --reads "data/raw/*_R{1,2}.fastq" --output "test"
```

### Other Parameters

#### Set custom parameters
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
