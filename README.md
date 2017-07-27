Overview
--------

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A50.25.1-brightgreen.svg)](https://www.nextflow.io/)

Influenza is an infectious disease caused by RNA viruses within the Orthomyxoviridae family. These viruses cause disease in a variety of animals including birds, pigs and humans. Because of their antigenically variable nature, these viruses are able to escape the innate and adaptive immune system to cause disease -- leading to the development of seasonal vaccine strains.

With the affordability, accessibility, and high-throughput nature of next generation sequencing instruments, individual labs are now able to sequence the genomes of a variety of organisms in bulk. In the influenza domain, this has led to a number of surveillance ([CDC](https://www.cdc.gov/flu/weekly/fluactivitysurv.htm), [Nextflu](https://academic.oup.com/bioinformatics/article/31/21/3546/194488/nextflu-real-time-tracking-of-seasonal-influenza), [USDA](https://www.aphis.usda.gov/aphis/ourfocus/animalhealth/animal-disease-information/avian-influenza-disease/ct_avian_influenza_disease)), prediction ([Goolge Flu Trends](http://people.sc.fsu.edu/~pbeerli/classes/ISC4931/ISC4931/SciComp/Entries/2013/2/18_Google_searches_and_influenza_files/detecting-influenza-epidemics.pdf), [Twitter Improves Influenza Forecasting](http://currents.plos.org/outbreaks/article/twitter-improves-influenza-forecasting/)), and sequence collection efforts ([Influenza Research Database](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5210613/), [Influenza Virus Database](https://www.ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi?go=genomeset)).

The purpose of this pipeline is to act as a preliminary analysis platform for the analysis of large amounts of influenza sequence data. The pipeline is built with [Nextflow](https://www.nextflow.io), a parallel computational workflow language. It also uses [Docker](https://www.docker.com), a software containerization platform that aids in the installation of the many open-source programs utilized throughout the pipeline. Functionally, the pipeline performs various quality control checks, host-dna removal, sequence assembly, complete genome annotation, genome coverage statistics, and easy-to-view html summary reports.

More Information
----------------
  - [Software Requirements](https://github.com/cdeanj/ai-assembly-pipeline/blob/master/docs/requirements.md)
  - [Process Overview](https://github.com/cdeanj/ai-assembly-pipeline/blob/master/docs/process.md)
  - [Installation](https://github.com/cdeanj/ai-assembly-pipeline/blob/master/docs/installation.md)
  - [Usage](https://github.com/cdeanj/ai-assembly-pipeline/blob/master/docs/usage.md)
  - [Configuration](https://github.com/cdeanj/ai-assembly-pipeline/blob/master/docs/configuration.md)
  - [Output](https://github.com/cdeanj/ai-assembly-pipeline/blob/master/docs/output.md)
  - [Dependencies](https://github.com/cdeanj/ai-assembly-pipeline/blob/master/docs/dependencies.md)
  - [Contact](https://github.com/cdeanj/ai-assembly-pipeline/blob/master/docs/contact.md)
  - [Acknowledgements](https://github.com/cdeanj/ai-assembly-pipeline/blob/master/docs/acknowledgements.md)
