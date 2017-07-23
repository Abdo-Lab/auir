Overview
--------

Influenza is an infectious disease caused by RNA viruses within the Orthomyxoviridae family. These viruses cause disease in a variety of animals including birds, pigs and humans. Because of their antigenically variable nature, these viruses are able to escape the innate and adaptive immune system to cause disease -- leading to the development of seasonal vaccine strains.

With the affordability, accessibility, and high-throughput nature of next generation sequencing instruments, individual labs are now able to sequence the genomes of a variety of organisms in bulk. In the influenza domain, this has led to a number of surveillance (Nextflu, CDC), prediction (Goolge Flu Trends, Twitter), and sequence collection efforts (Influenza Research Database, Influenza Virus Database).

The purpose of this pipeline is to act as a prelimnary analysis platform for the analysis of large amounts of influenza sequence data. The pipeline is built with Nextflow, a parallel computational workflow language. It also uses Docker, a software containerization platform that aids in the installation of the many open-source programs utilized throughout the pipeline. Functionally, the pipeline performs various quality control metrics, sequence assembly, complete genome annotations, genome coverage statistics, and easy-to-view html summary reports.

More Information
----------------

  - Installation
  - Usage
  - Configuration
  - Output
  - Dependencies
  - Acknowledgements
