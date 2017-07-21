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
 https://github.com/cdeanj/ai-assembly-pipeline
--------------------------------------------------------------------------------
 @Documentation
 https://github.com/cdeanj/ai-assembly-pipeline/blob/master/README.md
--------------------------------------------------------------------------------
@Licence
 https://github.com/cdeanj/ai-ssembly-pipeline/blob/master/LICENSE
--------------------------------------------------------------------------------
 Processes overview
 - Run AI pipeline
================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/

params.reads = "$baseDir/data/raw/*_R{1,2}_001.fastq"
params.adapters = "$baseDir/data/adapters/nextera.fa"
params.output = "./test"
params.threads = 1

adapters = file(params.adapters)
threads = params.threads

params.leading = 3
params.trailing = 3
params.slidingwindow = "4:15"
params.minlen = 36

leading = params.leading
trailing = params.trailing
slidingwindow = params.slidingwindow
minlen = params.minlen

Channel
    .fromFilePairs( params.reads, flat: true )
    .ifEmpty { exit 1, "Read pair files could not be found: ${params.reads}" }
    .into { read_pairs }

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
        set dataset_id, file("${dataset_id}.trimmomatic.stats.log") into (log)

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

    publishDir "${params.output}/SPAdes_Contigs", mode: "copy"

    input:
        set dataset_id, file(forward), file(reverse) from paired_fastq

    output:
        set dataset_id, file("${dataset_id}.contigs.fa") into (spades_contigs)

    """
    spades.py \
      -t ${threads} \
      --only-assembler \
      -1 ${forward} \
      -2 ${reverse} \
      -o output

    mv output/contigs.fasta .
    mv contigs.fasta ${dataset_id}.contigs.fa
    """
}

/*process RemoveMinContigs {
	publishDir "${params.output}/SPAdes_Contigs", mode: "copy"

	maxForks 8

	tag { dataset_id }

	input:
	set dataset_id, file(contigs) from spades_contigs

	output:
	set dataset_id, file("${dataset_id}_contigs.fa") into min_contigs

	"""
        #!/usr/bin/python       
        from Bio import SeqIO

        fasta_fp = "${contigs}"
        record_dict = list(SeqIO.parse(fasta_fp, "fasta"))

        target_records = []
        for record in record_dict:
                if len(record.seq) >= 200:
                        target_records.append(record)

        SeqIO.write(target_records, "${dataset_id}_contigs.fa", "fasta")
	"""
}*/

/*process BlastContigs {
	publishDir "${params.output}/Blast", mode: "copy"

	tag { dataset_id }

	input:
	set dataset_id, file(contigs) from min_contigs

	output:
	set dataset_id, file("${dataset_id}_annotated_contigs.fa") into annotated_contigs

	"""
	blastn -db nt -query ${contigs} -num_alignments 1 -outfmt "10 stitle" -num_threads ${threads} > ${dataset_id}_annotations
	cat ${dataset_id}_annotations | sed -e '/Influenza/s/^/>/' > ${dataset_id}_annotations.txt
	awk '/^>/ { getline <"${dataset_id}_annotations.txt" } 1 ' ${contigs} > ${dataset_id}_annotated_contigs.fa
	"""
}*/

/*process RemoveDuplicateAnnotations {
	publishDir "${params.output}/CleanedContigs", mode: "copy"

	input:
	set dataset_id, file(contigs) from annotated_contigs

	output:
	set dataset_id, file("${dataset_id}_clean_contigs.fa") into cleaned_contigs

	"""
	#!/usr/bin/python
	from Bio import SeqIO

	fasta_fp = "${contigs}"
	record_dict = list(SeqIO.parse(fasta_fp, "fasta"))

	mapper = {}

	for record in record_dict:
		if record.description in mapper:
			if len(record.seq) > mapper[record.description]:
				mapper[record.description] = len(record.seq)

		else:
			mapper[record.description] = len(record.seq)

	target = []

	for record in record_dict:
		if record.description in mapper and len(record.seq) == mapper[record.description]:
			target.append(record)

	SeqIO.write(target, "${dataset_id}_clean_contigs.fa", "fasta")
	"""
}*/



/*process QuastEvaluation {
	publishDir "${params.output}/QuastEvaluation", mode: "move"

	tag { dataset_id }

	input:
	set dataset_id, file(quast_contigs) from cleaned_contigs

	output:
	file("output/*") into quast_evaluation
	
	"""
	quast.py ${quast_contigs} --space-efficient --threads ${threads} -o output
	"""
}*/




/*process AlignReadsToContigs {
	publishDir "${params.output}/Alignment", mode: "copy"

	input:
	set dataset_id, file()

	output:
	set dataset_id, file("${dataset_id}_alignment.sam")

	"""
	bwa index ${contigs}
	bwa mem ${contigs} ${forward} ${reverse} > ${dataset_id}_alignment.sam
	"""
}*/


/*process AlignToReference {
	publishDir "${params.output}/Alignment", mode: "copy"

	tag { dataset_id }

	input:
	set dataset_id, file(contigs) from min_contigs
	file genome
	file index from genome_index.first()

	output:
	set dataset_id, file("${dataset_id}_alignment.sam") into quast_contigs

	"""
	bwa mem -t ${threads} ${genome} ${contigs} > ${dataset_id}_alignment.sam
	"""
}*/

