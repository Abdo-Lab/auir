Installation
------------

To run the pipeline, you will need Nextflow and the appropriate Docker containers to run each process within the pipeline.

Install Nextflow
----------------
```
$ curl -fsSL get.nextflow.io | bash
$ ./nextflow
$ mv nextflow /usr/local/bin
```

Pull Docker Containers
----------------------
This can take anywhere between 10 and 15 minutes depending on your connection speed.
```
$ docker pull cdeanj/ai
```
