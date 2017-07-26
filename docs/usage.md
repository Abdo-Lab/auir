Usage
-----

### Display Help Message
```
$ nextflow run auir.nf --help
```

```
Usage: 
  ./nextflow run auir.nf --reads <reads> --adapters <adapters> --threads <threads> --output <output>
Where: 
     <reads> is the directory location of your FASTQ formatted read pairs
  <adapters> is the location of your FASTA formatted adapter sequences
   <threads> is the number of threads to use for each process
    <output> is the output directory to write intermediate outputs to
     <index> is the directory location of your host genome index files
      <host> is the location of your FASTA formatted host genome
      <help> displays the help message
Contact: 
  Christopher Dean <cdean11@colostate.edu>
Version: 
  Auir 1.0
```
-----
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
-----
### File Outputs

```
$ nextflow run auir.nf --reads "data/raw/*_R{1,2}.fastq" --output "test"
```
-----
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
-----



