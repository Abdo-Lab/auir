#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
================================================================================
=                           A I - c o n t a i n e r s                          =
================================================================================
@Author
Christopher Dean <cdean11@colostate.edu>
--------------------------------------------------------------------------------
 @Homepage
 https://github.com/cdeanj/ai-containers
--------------------------------------------------------------------------------
 @Documentation
 https://github.com/cdeanj/ai-containers/blob/master/README.md
--------------------------------------------------------------------------------
@Licence
 https://github.com/cdeanj/ai-containers/blob/master/LICENSE
--------------------------------------------------------------------------------
 Processes overview
 - Run AI pipeline
================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/

params.reads = "$baseDir/cov_data/*_R{1,2}_001.fq"
//params.reads = "/home/chris_dean/projects/influenza/scripts/14999-1_S5_L001_R1_001/*_R{1,2}_001.fq"
params.adapters = "$baseDir/adapters/adapters.fa"
params.output = "./TEST"
params.host = ""
params.genome = "$baseDir/data/genome/H3N8.fasta"
params.host_index = ""
params.genome_index = ""
params.threads = 1;

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


/*
 * Validate input parameters
 */
Channel
        .fromFilePairs( params.reads, flat: true )
	.ifEmpty { exit 1, "Read pairs could not be found: ${params.reads}" }
        .into { read_pairs }

if( params.host ) {
	host = file(params.host)
	if( !host.exists() ) exit 1, "Host genome could not be found ${params.host}"
}

if( params.genome) {
	genome = file(params.genome)
	if( !genome.exists() ) exit 1, "Reference genome could not be found ${params.genome}"
}

if( params.host_index ) {
	host_index = Channel.fromPath( params.host_index )
	if( !host_index.exists() ) exit 1, "Host genome index files could not be found ${params.host_index}"
}

if( params.genome_index ) {
	genome_index = Channel.fromPath( params.genome_index )
	if( !genome_index.exists() ) exit 1, "Reference genome index files could not be found ${params.genome_index}"
}

/*
 * Remove low quality base pairs and adapter sequences with Trimmomatic
 */
process QC {
	container 'chrisd/virus'

	publishDir "${params.output}/Preprocessing", mode: "copy"

	errorStrategy 'ignore'
	
	maxForks 8

	tag { dataset_id }

	input:
	set dataset_id, file(forward), file(reverse) from read_pairs

	output:
	set dataset_id, file("${dataset_id}_1P.fastq"), file("${dataset_id}_2P.fastq") into (trimmed_read_pairs, bwa_read_pairs)

	"""
	java -jar ${TRIMMOMATIC}/trimmomatic-0.36.jar PE \
	-threads ${threads} \
	$forward $reverse -baseout ${dataset_id} \
	ILLUMINACLIP:Trimmomatic-0.36/adapters/${adapters}:2:30:10:3:TRUE \
	LEADING:${leading} \
	TRAILING:${trailing} \
	SLIDINGWINDOW:${slidingwindow} \
	MINLEN:${minlen}
	mv ${dataset_id}_1P ${dataset_id}_1P.fastq
        mv ${dataset_id}_2P ${dataset_id}_2P.fastq
	"""
}


/*
 * Generate genome assemblies with SPAdes
 */
process SPAdesAssembly {
	container 'chrisd/virus'

	publishDir "${params.output}/SPAdes_Contigs", mode: "copy"

	errorStrategy 'ignore'

	maxForks 8

	tag { dataset_id }
	
	input:
	set dataset_id, file(forward), file(reverse) from trimmed_read_pairs

	output:
	set dataset_id, file("${dataset_id}_contigs.fa") into spades_contigs

	"""
	spades.py -t ${threads} --pe1-1 ${forward} --pe1-2 ${reverse} -o output
	mv output/contigs.fasta .
	mv contigs.fasta ${dataset_id}_contigs.fa
	"""
}

process RemoveMinContigs {
	container 'chrisd/virus'

	errorStrategy 'ignore'

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
}

process BlastContigs {
	publishDir "${params.output}/Blast", mode: "copy"

	errorStrategy 'ignore'

	maxForks 8

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
}

process RemoveDuplicateAnnotations {
	publishDir "${params.output}/CleanedContigs", mode: "copy"

	errorStrategy 'ignore'

	maxForks 8

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
}



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
	set dataset_id, file(

	output:
	set dataset_id, file("${dataset_id}_alignment.sam")

	"""
	bwa index ${contigs}
	bwa mem ${contigs} ${forward} ${reverse} > ${dataset_id}_alignment.sam
	"""
}*/


/*process AlignToReference {
	errorStrategy 'ignore'

	publishDir "${params.output}/Alignment", mode: "copy"

	maxForks 8

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

