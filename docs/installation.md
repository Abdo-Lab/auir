Installation
------------

To run the pipeline, you will need Nextflow and the appropriate Docker containers to run each process within the pipeline.

### Install Nextflow
```
$ curl -fsSL get.nextflow.io | bash
$ ./nextflow
$ mv nextflow /usr/local/bin
```

### Pull Docker Containers
```
$ docker pull cdeanj/auir -a 
```
This command will use Docker to download each of the tools required to run the pipeline. The Dockerfiles for each of downloaded tools can be found on [Dockerhub](https://hub.docker.com/r/cdeanj/auir/). This can take between 10 and 15 minutes depending on your connection speed.

### Clone Github Repository
```
$ git clone https://github.com/cdeanj/auir.git
$ cd auir
```
