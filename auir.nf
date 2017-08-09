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


if( params.help ) {
    return help()
}
if( !nextflow.version.matches('0.25+') ) {
    return nextflow_version_error()
}
if( params.index ) { 
    index = Channel.fromPath(params.index).toSortedList() 
    if( !index.exists() ) return index_error(index)
}
if( params.host ) {     
    host = file(params.host)
    if( !host.exists() ) return host_error(host)
}
if( params.adapters ) {     
    adapters = file(params.adapters) 
    if( !adapters.exists() ) return adapter_error(adapters)
}
if( params.fqc_adapters ) {
    fqc_adapters = file(params.fqc_adapters)                             
    if( !fqc_adapters.exists() ) return fastqc_error(fqc_adapters)
}

threads = params.threads

leading = params.leading
trailing = params.trailing
slidingwindow = params.slidingwindow
minlen = params.minlen

min_contig = params.min_contig

Channel
    .fromFilePairs( params.reads, flat: true )
    .ifEmpty { return fastq_error(params.reads) }
    .into { read_pairs; fastqc_pairs; alignment_pairs }

process RunFastQC {
    tag { dataset_id }

    publishDir "${params.output}/RunFastQC", mode: 'copy'

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

process RunQC {
    tag { dataset_id }

    publishDir "${params.output}/RunQC", mode: 'copy',
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

if( !params.index ) {
    process BuildHostIndex {
        tag { host.baseName }

        publishDir "${params.output}/BuildHostIndex", mode: "copy"

        input:
            file(host)

        output:
            file '*' into index

        """
        bwa index ${host}
        """
    }
}

process AlignReadsToHost {
    tag { host.baseName }
        
    publishDir "${params.output}/AlignReadsToHost", mode: "copy"
        
    input:
        set dataset_id, file(forward), file(reverse) from paired_fastq
        file idx from index.first()
        file host
            
    output:
        set dataset_id, file("${dataset_id}.host.sam") into host_sam
            
    """ 
    bwa mem ${host} ${forward} ${reverse} -t ${threads} > ${dataset_id}.host.sam
    """ 
}

process RemoveHostDNA {
    tag { dataset_id }

    publishDir "${params.output}/RemoveHostDNA", mode: "copy"

    input:
        set dataset_id, file(sam) from host_sam

    output:
        set dataset_id, file("${dataset_id}.host.sorted.removed.bam") into host_bam

    """
    samtools view -bS ${sam} | samtools sort -@ ${threads} -o ${dataset_id}.host.sorted.bam
    samtools view -h -f 4 -b ${dataset_id}.host.sorted.bam -o ${dataset_id}.host.sorted.removed.bam
    """
}

process BAMToFASTQ {
    tag { dataset_id }

    publishDir "${params.output}/BAMToFASTQ", mode: "copy"

    input:
        set dataset_id, file(bam) from host_bam

    output:
        set dataset_id, file("${dataset_id}.non.host.R1.fastq"), file("${dataset_id}.non.host.R2.fastq") into non_host_fastq

    """
    bedtools  \
       bamtofastq \
      -i ${bam} \
      -fq ${dataset_id}.non.host.R1.fastq \
      -fq2 ${dataset_id}.non.host.R2.fastq
    """
}

process RunSPAdes {
    tag { dataset_id }

    publishDir "${params.output}/RunSPAdes", mode: "copy"

    input:
        set dataset_id, file(forward), file(reverse) from non_host_fastq

    output:
        set dataset_id, file("${dataset_id}.contigs.fa") into (spades_contigs)

    script:
    """
    spades.py \
      -t ${threads} \
      --only-assembler \
      --cov-cutoff auto \
      -1 ${forward} \
      -2 ${reverse} \
      -o output

    mv output/contigs.fasta .
    # Remove line breaks from FASTA files
    # Taken from Kent at https://stackoverflow.com/questions/15857088/remove-line-breaks-in-a-fasta-file
    awk '/^>/ { print (NR==1 ? "" : RS) \$0; next } { printf "%s", \$0 } END { printf RS }' contigs.fasta > ${dataset_id}.contigs.temp.fa
    python $baseDir/bin/min_contigs.py -r ${dataset_id}.contigs.temp.fa -m ${min_contig} -o ${dataset_id}.contigs.min.removed.fa
    python $baseDir/bin/dupe_contigs.py -r ${dataset_id}.contigs.min.removed.fa -o ${dataset_id}.contigs.fa
    """
}

process RunBlast {
    tag { dataset_id }

    publishDir "${params.output}/RunBlast", mode: 'copy'

    input:
        set dataset_id, file(contigs) from spades_contigs

    output:
        set dataset_id, file("${dataset_id}.contigs.annotated.fa") into (annotated_spades_contigs, annotated_spades_contigs2, annotated_spades_contigs3)

    """
    blastn -db InfluenzaDB -query ${contigs} -max_hsps 1 -max_target_seqs 1 -outfmt "10 stitle" -num_threads ${threads} -task megablast > ${dataset_id}.contig.description.tmp
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

process AlignReadsToContigs {
    tag { dataset_id }

    publishDir "${params.output}/AlignReadsToContigs", mode: "copy"

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

    publishDir "${params.output}/SAMToBAM", mode: "copy"

    input:
        set dataset_id, file(sam) from sams

    output:
        set dataset_id, file("${dataset_id}.alignment.bam") into (bams)

    """
    samtools view -bS ${sam} | samtools sort -@ ${threads} -o ${dataset_id}.alignment.bam
    """
}

process RemovePCRDuplicates {
    tag { dataset_id }

    publishDir "${params.output}/RemovePCRDuplicates", mode: "copy", pattern: "*.bam"

    input:
        set dataset_id, file(bam) from bams

    output:
        set dataset_id, file("${dataset_id}.alignment.rmdup.bam"), file("${dataset_id}.alignment.rmdup.bam.bai") into (pcr_rmdup_bams, bedtools_bams)

    """
    samtools rmdup ${bam} ${dataset_id}.alignment.rmdup.bam
    samtools index ${dataset_id}.alignment.rmdup.bam
    """
}

freebayes_tuple = pcr_rmdup_bams.combine(annotated_spades_contigs3, by: 0)

process RunFreebayes {
    tag { dataset_id }

    publishDir "${params.output}/RunFreebayes", mode: "copy"

    input:
        set dataset_id, file(bam), file(bai), file(contigs) from freebayes_tuple

    output:
        set dataset_id, file("${dataset_id}.vcf") into vcfs

    """
    freebayes -f ${contigs} -p 2 -C 5 ${bam} > ${dataset_id}.vcf
    """
}

process CalculateCoverage {
    tag { dataset_id }

    publishDir "${params.output}/CalculateCoverage", mode: "copy"

    input:
        set dataset_id, file(bam) from bedtools_bams

    output:
        file("${dataset_id}.depth.coordinates") into (coordinates)

    """
    bedtools genomecov -d -ibam ${bam} > ${dataset_id}.depth.coordinates
    sed -i 's/^/${dataset_id}\t/' ${dataset_id}.depth.coordinates
    """
}

process AggregateCoverageCounts {
    publishDir "${params.output}/AggregateCoverageCounts", mode: "copy"

    input:
        file(counts) from coordinates.toSortedList()

    output:
        file("aggregated_counts.tsv") into aggregated_counts

    """
    cat ${counts} > aggregated_counts.tsv
    """
}

process RunQuast {
    input:
        file(annotated_contigs) from annotated_assemblies

    output:
        file("report.tsv") into (quast_logs)

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

process RunMultiQC {
    publishDir "${params.output}/RunMultiQC", mode: 'copy'

    input:
        file('*') from multiQCReports

    output:
        file("*multiqc_report.html") into multiQCReport

    """
    multiqc -f -v .
    """
}

def nextflow_version_error() {
    println ""
    println "This workflow requires Nextflow version 0.25 or greater -- You are running version $nextflow.version"
    println "Run ./nextflow self-update to update Nextflow to the latest available version."
    println ""
    return 1
}

def adapter_error(def input) {
    println ""
    println "[params.adapters] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def fastqc_error(def input) {
    println ""
    println "[params.fqc_adapters] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def fastq_error(def input) {
    println ""
    println "[params.reads] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def host_error(def input) {
    println ""
    println "[params.host] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def index_error(def input) {
    println ""
    println "[params.index] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def help() {
    println ""
    println "Program: Auir"
    println "Version: $workflow.repository - $workflow.revision [$workflow.commitId]"
    println "Documentation: https://github.com/cdeanj/auir"
    println "Contact: Christopher Dean <cdean11@colostate.edu>"
    println ""
    println "Usage:    nextflow run auir.nf [options]"
    println ""
    println "Input/output options:"
    println ""
    println "    --reads         STR      path to FASTQ formatted input sequences"
    println "    --adapters      STR      path to FASTA formatted adapter sequences"
    println "    --host          STR      path to FASTA formatted host genome"
    println "    --index         STR      path to BWA generated index files"
    println "    --output        STR      directory to write process outputs to"
    println ""
    println "Trimming options:"
    println ""
    println "    --leading       INT      cut bases off the start of a read, if below a threshold quality"
    println "    --minlen        INT      drop the read if it is below a specified length"
    println "    --slidingwindow INT      perform sw trimming, cutting once the average quality within the window falls below a threshold"
    println "    --trailing      INT      cut bases off the end of a read, if below a threshold quality"
    println ""
    println "Algorithm options:"
    println ""
    println "    --threads       INT      number of threads to use for each process"
    println "    --ploidy        INT      genome copy number"
    println "    --min-alt-count INT      requires this number of observations supporting an alternate allele"
    println ""
    println "Help options:"
    println ""
    println "    --help                   display this message"
    println ""
    return 1
}
