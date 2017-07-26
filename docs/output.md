Output
------

### Output Modules

Output files are written to the user-defined output directory (default: test/) Below we discuss the types of output files produced from each pipeline module and where they can be found on your system.

#### Module FastQC

FastQC produces html reports for raw sequence data in FASTQ format. These reports include content related to:
  - Sequence content
  - Sequence quality
  - GC content
  - Sequence length distribution
  - Adapter contamination
  
    
#### Module Trimmomatic

Trimmomatic removes low quality base pairs, adapter sequences and other sequence fragments. This module takes in as input a collection of FASTQ read pairs and a FASTA formatted adapter sequence file. By default, Nextera adapter sequences are used. 

As output, this module produces four outputs:
  - sample.1P.fastq
  - sample.2P.fastq
  - sample.1U.fastq
  - sample.2U.fastq

#### Module BWA

BWA is a sequence aligner that takes short fragments of DNA (or reads) and aligns them to a reference genome. It produces output in the standard SAM format and can be used by other tools to extract information related to the structure and quality of an alignment. This module serves two purposes:
  - Host-DNA removal
  - Per-base sequence coverage
  
In the host removal step, reads are aligned to the user-input host genome. This produces a SAM formatted alignment file that is used by Samtools to remove reads that aligned to the host.

To calculate per-base sequence coverage, reads 

#### Module SPAdes

SPAdes is a de-bruijn graph based genome assembler that takes short reads produced from, for example, Illumina instruments and breaks them into smaller pieces called (kmers) to produce a graph that can be traversed with the goal of reproducing the original genome. As input, this module takes in the trimmed FASTQ sequence files from Trimmomatic.

As output, this module procudes one output:
  - sample.contigs.fa

#### Module Blastn

#### Module Bedtools2

#### Module QUAST

QUAST is a tool for evaluating different features of a genome assembly and is a great way to guage the quality of an assembly.

#### Module MultiQC

MultiQC is a tool that can take outputs produced from open-source bioinformatics tools (log files, genome assemblies, or alignment files) and turn them into easy-to-view html reports.
