Installation
------------

To run the pipeline, you will need Nextflow and the appropriate Docker containers to run each process within the pipeline.

### Step 1 -- Download Nextflow
```
$ curl -fsSL get.nextflow.io | bash
$ ./nextflow
$ mv nextflow /usr/local/bin
```

### Step 2 -- Download Docker Containers

This command will download each of the required tools to run the pipeline. The Dockerfiles for each tool can be found on [Dockerhub](https://hub.docker.com/r/cdeanj/auir/).
```
$ docker pull cdeanj/auir -a 
```

When the download is complete, you should have nine Docker images.
```
$ docker images
REPOSITORY          TAG                 IMAGE ID            CREATED             SIZE
cdeanj/auir         bedtools2           7db588a1bb8e        About an hour ago   330.3 MB
cdeanj/auir         fastqc              27b4b5000102        About an hour ago   749.2 MB
cdeanj/auir         makeblastdb         3c196b1d1e92        About an hour ago   1.609 GB
cdeanj/auir         ivd                 6d5a5c9cb77c        About an hour ago   1.059 GB
cdeanj/auir         quast               a798c88e7b4a        About an hour ago   1.103 GB
cdeanj/auir         trimmomatic         010deabfa8e7        About an hour ago   104.9 MB
cdeanj/auir         spades              d5b62b86822f        About an hour ago   164.9 MB
cdeanj/auir         samtools            09ba9b5ec0b5        About an hour ago   313.1 MB
cdeanj/auir         bwa                 9d052c06d12b        About an hour ago   262.2 MB
```

### Step 3 -- Download Github Repository
```
$ git clone https://github.com/cdeanj/auir.git
$ cd auir
```
