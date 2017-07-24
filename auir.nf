#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
================================================================================
=                     A I | A S S E M B L Y | P I P E L I N E                  =
================================================================================
@Author
Christopher Dean <cdean11@colostate.edu>
--------------------------------------------------------------------------------
 @Homepage
 https://github.com/cdeanj/auir
--------------------------------------------------------------------------------
 @Documentation
 https://github.com/cdeanj/auir/blob/master/README.md
--------------------------------------------------------------------------------
 @Licence
 https://github.com/cdeanj/auir/blob/master/LICENSE
--------------------------------------------------------------------------------
 Processes overview
 - Run auir influenza pipeline
================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/

threads = params.threads
adapters = file(params.adapters)
fqc_adapters = file(params.fqc_adapters)
leading = params.leading
trailing = params.trailing
slidingwindow = params.slidingwindow
minlen = params.minlen

Channel
    .fromFilePairs( params.reads, flat: true )
    .ifEmpty { exit 1, "Read pair files could not be found: ${params.reads}" }
    .into { read_pairs; fastqc_pairs; alignment_pairs }

process FastQC {
    tag { dataset_id }

    publishDir "${params.output}/FastQC", mode: 'copy'

    input:
        set dataset_id, file(forward), file(reverse) from fastqc_pairs

    output:
        set dataset_id, file("*_fastqc.zip") into (fastqc_logs)

    """
    mkdir output
    fastqc -f fastq ${forward} ${reverse} -t ${threads} -o output
    mv output/*.zip .
    """
}

process Trimmomatic {
    tag { dataset_id }

    publishDir "${params.output}/Trimmomatic", mode: 'copy',
        saveAs: { filename ->
            if(filename.indexOf("P.fastq") > 0) "Paired/$filename"
            else if(filename.indexOf("U.fastq") > 0) "Unpaired/$filename"
            else if(filename.indexOf(".log") > 0) "Log/$filename"
            else {}
        }
	
    input:
        set dataset_id, file(forward), file(reverse) from read_pairs

    output:
        set dataset_id, file("${dataset_id}.1P.fastq"), file("${dataset_id}.2P.fastq") into (paired_fastq)
        set dataset_id, file("${dataset_id}.1U.fastq"), file("${dataset_id}.2U.fastq") into (unpaired_fastq)
        set dataset_id, file("${dataset_id}.trimmomatic.stats.log") into (trimmomatic_logs)

    """
    java -jar ${TRIMMOMATIC}/trimmomatic-0.36.jar \
      PE \
      -threads ${threads} \
      $forward $reverse -baseout ${dataset_id} \
      ILLUMINACLIP:${adapters}:2:30:10:3:TRUE \
      LEADING:${leading} \
      TRAILING:${trailing} \
      SLIDINGWINDOW:${slidingwindow} \
      MINLEN:${minlen} \
      2> ${dataset_id}.trimmomatic.stats.log

    mv ${dataset_id}_1P ${dataset_id}.1P.fastq
    mv ${dataset_id}_2P ${dataset_id}.2P.fastq
    mv ${dataset_id}_1U ${dataset_id}.1U.fastq
    mv ${dataset_id}_2U ${dataset_id}.2U.fastq
    """
}

process SPAdes {
    tag { dataset_id }

    publishDir "${params.output}/SPAdes", mode: "copy"

    input:
        set dataset_id, file(forward), file(reverse) from paired_fastq

    output:
        set dataset_id, file("${dataset_id}.contigs.fa") into (spades_contigs)

    """
    spades.py \
      -t ${threads} \
      --only-assembler \
      --cov-cutoff auto \
      -1 ${forward} \
      -2 ${reverse} \
      -o output

    mv output/contigs.fasta .
    mv contigs.fasta ${dataset_id}.contigs.fa
    """
}

process Blastn {
    tag { dataset_id }

    publishDir "${params.output}/Blastn", mode: 'copy'

    input:
        set dataset_id, file(contigs) from spades_contigs

    output:
        set dataset_id, file("${dataset_id}.contigs.annotated.fa") into (annotated_spades_contigs, annotated_spades_contigs2)

    """
    blastn -db InfluenzaDB -query ${contigs} -max_hsps 1 -max_target_seqs 1 -outfmt "10 stitle" -num_threads ${threads} > ${dataset_id}.contig.description.tmp
    cat ${dataset_id}.contig.description.tmp | sed -e '/Influenza/s/^/>/' > ${dataset_id}.contig.description.txt
    sed -i 's/ /_/g' ${dataset_id}.contig.description.txt
    awk '/^>/ { getline <"${dataset_id}.contig.description.txt" } 1 ' ${contigs} > ${dataset_id}.contigs.annotated.fa
    """
}

annotated_assemblies = Channel.empty()
    .mix(annotated_spades_contigs)
    .flatten()
    .toList()

reads_and_contigs = alignment_pairs.combine(annotated_spades_contigs2, by: 0)

process BWA {
    tag { dataset_id }

    publishDir "${params.output}/BWA", mode: "copy"

    input:
        set dataset_id, file(forward), file(reverse), file(contigs) from reads_and_contigs

    output:
        set dataset_id, file("${dataset_id}.alignment.sam") into (sams)

    """
    bwa index ${contigs}
    bwa mem ${contigs} ${forward} ${reverse} -t ${threads} > ${dataset_id}.alignment.sam
    """
}

process SAMToBAM {
    tag { dataset_id }

    publishDir "${params.output}/BWA", mode: "copy"

    input:
        set dataset_id, file(sam) from sams

    output:
        set dataset_id, file("${dataset_id}.alignment.bam") into (bams)

    """
    samtools view -bS ${sam} | samtools sort -@ ${threads} -o ${dataset_id}.alignment.bam
    """
}

process Bedtools {
    tag { dataset_id }

    publishDir "${params.output}/Bedtools", mode: "copy"

    input:
        set dataset_id, file(bam) from bams

    output:
        file("${dataset_id}.depth.coordinates") into (coordinates)

    """
    bedtools genomecov -d -ibam ${bam} > ${dataset_id}.depth.coordinates
    sed -i 's/^/${dataset_id}\t/' ${dataset_id}.depth.coordinates
    """
}

process AggregateCounts {
    publishDir "${params.output}/Bedtools", mode: "copy"

    input:
        file(counts) from coordinates.toSortedList()

    output:
        file("aggregated_counts.tsv") into res

    """
    cat ${counts} > aggregated_counts.tsv
    """
}



/*process QUAST {
    publishDir "${params.output}/QUAST", mode: 'copy'

    input:
        file(annotated_contigs) from annotated_assemblies

    output:
        set file("report.tsv") into (quast_logs)

    """
    quast.py \
      ${annotated_contigs} \
      --no-plots \
      --no-html \
      --no-icarus \
      --no-snps \
      --no-sv \
      --est-ref-size 13600 \
      -t ${threads} \
      -o output

    mv output/report.tsv .
    """
}

multiQCReports = Channel.empty()
    .mix(
        trimmomatic_logs,
        quast_logs,
        fastqc_logs
    )
    .flatten().toList()

process MultiQC {
    publishDir "${params.output}/MultiQC", mode: 'copy'

    input:
        file('*') from multiQCReports

    output:
        set file("*multiqc_report.html") into multiQCReport

    """
    multiqc -f -v .
    """
}*/
