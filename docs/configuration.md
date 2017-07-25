Configuration
-------------

The pipeline is configured with a configuration file that can be used to modify input parameters for a variety of processes. To modify these parameters, open the nextflow.config file and update the desired parameters.

### Input Sequences
```
/* Location of forward and reverse read pairs */
reads = "data/raw/*_R{1,2}_001.fastq"
```

### Index Files
```
/* Location of host genome index files */
index = "data/genome/index/*"
```

### Host Genome
```
/* Location of host genome file */
host = "data/genome/host/gallus.fa"
```

### Adapter Sequences
```
/* Location of adapter sequences */
adapters = "data/adapters/nextera.fa"
```

### Output Directory
```
/* Output directory */
output = "./test"
```

### Number of Threads
```
/* Number of threads */
threads = 16
```

### Trimmomatic Parameters
```
/* Trimmomatic trimming parameters */
leading = 3
trailing = 3
slidingwindow = "4:15"
minlen = 36
```
