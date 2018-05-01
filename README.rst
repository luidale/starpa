starpa
======

.. image:: https://img.shields.io/pypi/v/starpa.svg
    :target: https://pypi.python.org/pypi/starpa
    :alt: Latest PyPI version

.. image:: https://travis-ci.org/luidale/starpa.png
   :target: https://travis-ci.org/luidale/starpa
   :alt: Latest Travis CI build status

**Stable RNA processing product analyzer**

Tool to predict, quantify and characterize stable RNA processing products
from RNA-seq data.

Overview
--------
Starpa workflow is divided into multiple consecutive tasks which can be executed separately, 
as a freely chosen successive subsets or all tasks at once in sequential order.
This adds flexibility to the tool to use as an input RNA-seq data in various state of processing.
For example Starpa can handle raw data in FastQ format, but also trimmed reads (FastQ format)
or aligned reads in SAM format.

Both paired-end (PE) and single-end (SE) sequencing reads are accepted as an input.

In addition, the tool is highly configurable and can handle multiple libraries in parallel manner (multiprocessing).

**Tasks are following:**

- *trim*

Cutadapt is used to trim low quality 3' end of the reads followed by adapter removal from 3' end 
of the reads. 

In case of SE, the reads where 3' adapter was not trimmed are excluded. 
This ensures that 3' end of the read is stable RNA processing products is estimated with higher 
confidence.

- *align*

Bowtie2 is used to align reads to the genome. All matches to the genome are recorded.

- *sam_sort*

From aligned reads the unmapped and discordantly mapped reads are discarded. In addition, only the reads belonging to 
best stratum (class of alignment score) are retained while alignments with lower alignments score 
are excluded.

- *pseudoSE*

Alignments with too many mismatches and reads with too many genomic alignments are discarded.
All other reads get NH tag (if not present) describing the number of reported alignments. 
Sequence and quality fields of secondary alignments are filled with sequence and quality data.
In the end the PE reads are converted to pseudo SE reads to ease subsequent analysis steps. 

- *identify*

Flaimapper2 is used to predict stable RNA processing products. To ensure prediction of all
processing products which share start or end positions, the reads are fractionated according 
to their length. Subsequently, Flaimmper2 is run on each fraction of reads separately and 
the predicted processing products are filtered by the read count (estimation by 
Flaimapper-2) exceeding threshold set. The filtered predicted processing products are quantified 
more precisely via bedtools intersect.

- *cluster*

Quantified processing products are filtered once again by the read counts (bedtools intersect)
exceeding threshold and by relative coverage (average coverage of reads assigned to processing products 
divided by average coverage of all reads aligned to the positions of processing products).
Next, the processing products from all libraries analysed are combined (identifying unique species) 
and clustered.

Clustering is two step process:

a) clustering by overlap.

As the prediction of processing products by Flaimapper-2 is probabilistic, the predicted ends 
of the processing products in different libraries might slightly vary, as also the true ends. 
Therefore, the predicted processing products which do largely overlap and have some bases 
(adjustable) not overlapping are clustered and representative processing products for clusters 
are selected.

b) clustering by sequence

As a majority of genomes contain repeating regions (repeat regions, rRNA operons, some tRNA genes etc)
reads can be mapped to multiple positions resulting multiple processing products consisting 
from the same or similar set of reads.
To reduce the number of identical processing products they are clustered by sequence identity 
via CDI-HIT-EST. Still the genomic matches of particular reads can be in genomic regions with different surrounding
sequence/context (eg. different genes) therefore clustering solely based on sequence identity can result 
loss of information.
To avoid it the predicted processing products which cluster by sequence identity has to be supported by the 
clustering (again via CDI-HIT-EST) of the contigs they overlap with and representative processing product for the 
clusters are selected.

In addition, the contigs are identified and wig formatted files (containing coverage data of 
individual libraries) are created.

- *quantify*

Representative processing products will be quantified using bedtools intersect in every library.
Additional characteristics will be gathered (relative coverage, coverage at single position level, 
consensus sequence, quality of consensus sequence, genomic sequence, uniqueness). Quantification data
is also converted to read per million of mapped reads (RPM), RPM of biotype and RPM of biotype groups.

Installation
------------
::

 pip install --user starpa


Requirements
^^^^^^^^^^^^
Starpa is depending on following tools which have to be installed in your system:

`Python3.4+ <https://www.python.org/>`_,
`cutadapt <https://github.com/marcelm/cutadapt>`_,
`bowtie2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_,
`samtools <http://www.htslib.org/doc/samtools.html>`_,
`Flaimapper-2 <hhttps://github.com/yhoogstrate/flaimapper>`_,
`bedtools <http://bedtools.readthedocs.io/en/latest/#>`_,
`CDI-HIT-EST <http://weizhongli-lab.org/cd-hit/>`_,
`featureCounts <http://bioinf.wehi.edu.au/featureCounts/>`_ (Release 1.6.1+).

Python3 requires following packages which will be installed (if missing) during 
the installation of starpa:

pyfaidx, docopt, schema

Compatibility
-------------
**OS:**

Starpa is compatible with UNIX like operating systems.

**Input:**

1) Colorspace reads are not supported.

2) Both paired-end (PE) and single-end (SE) reads are supported.

Usage
-----
Usage of starpa is as follows::

 Usage:
     starpa      [-hv]
     starpa -s <start_task> -e <end_task> -c <parameter_file> -i <input> 
     -o <output>

 Arguments:

     <start_task>        task to start with
     <end_task>          tast to end with
     <config_file>       configuration file
     <input>             input folder
     <output>            output folder
 Options:
     -v, --version
     -h, --help
     -s <start_task>, --start=<start_task>
     -e <end_task>, --end=<end_task>
     -c <config_file>, --config=<config_file>
     -i <input_folder>, --input=<input_folder>
     -o <output_folder>, --output=<output_folder>

|

**Tasks**

Starpa work-flow is divided into multiple consecutive tasks which can be executed:

- separately
- as a freely chosen successive subsets 
- all at once in sequential order

Tasks in sequential order:

	trim, align, sam_sort, pseudoSE, identify, cluster, quantify

**Configuration file**

`Configuration file <https://raw.githubusercontent.com/luidale/starpa/master/src/starpa/data/config.txt>`_ 
is used to set various parameters which allow to adjust the 
performance of the work-flow according to the user needs and input data.
The description of each parameter is given in the file itself.

Configuration file states also the location of following files:

adapter files - adapter sequencies in fasta format

genome file - genome sequence in fasta format

annotation file - in GFF or GFF3 format.

`"flaimapper parameter file" <https://raw.githubusercontent.com/luidale/starpa/master/src/starpa/data/flaimapper_parameters/parameters.dev-2-100-2.txt>`_  -
described in more deteil `here <https://github.com/yhoogstrate/flaimapper#the---parameters-argument>`_. Given Flaimapper-2 parameters file is adjusted to be suitable to predict processing products with rather defined ends.

`"library_file" <https://raw.githubusercontent.com/luidale/starpa/master/src/starpa/data/libraries.txt>`_ - 
describing libraries to be analysed.

"library_file" is a tabular file containing:
 1) the name of the libraries

 2) conditions they are derived from and 

 3) identifier of replicate 

(note that all three columns are separated by tab)

::

 #Library number	Sample	Replicate
 library1	LB OD 0.4	I
 library2	LB OD 0.4	II

| 

`Configuration file <https://raw.githubusercontent.com/luidale/starpa/master/src/starpa/data/config.txt>`_,
`"flaimapper parameter file" <https://raw.githubusercontent.com/luidale/starpa/master/src/starpa/data/flaimapper_parameters/parameters.dev-2-100-2.txt>`_ and
`"library_file" <https://raw.githubusercontent.com/luidale/starpa/master/src/starpa/data/libraries.txt>`_ are available in:

::

 src/starpa/data

|


**Input folder**

While running a single or multiple tasks, the input folder has to contain specific data 
required for the first task. 
For the following task the preceding tasks will prepare proper data.

Each task has different requirements for the input data:

- *trim*

| Sequencing data in `FastQ format <https://en.wikipedia.org/wiki/FASTQ_format>`_.
| Can be in PE or SE format which has to be indicated in 
 `configuration file <https://raw.githubusercontent.com/luidale/starpa/master/src/starpa/data/config.txt>`_ .
| FastQ files can be compressed as ".gz", ".bz2" or ".xz".


- *align*

| Trimmed and cleaned reads in `FastQ format <https://en.wikipedia.org/wiki/FASTQ_format>`_.
| Can be in PE or SE format which has to be indicated in 
 `configuration file <https://raw.githubusercontent.com/luidale/starpa/master/src/starpa/data/config.txt>`_ .
| FastQ files can be compressed as ".gz" (requires bowtie2.3.1+)


- *sam_sort*

| Aligned reads in SAM format. 
| Can be in PE or SE format which has to be indicated in 
 `configuration file <https://raw.githubusercontent.com/luidale/starpa/master/src/starpa/data/config.txt>`_ .

| BAM format is not currently supported.


- *pseudoSE*

| Aligned reads in SAM format. 
| Can be in PE or SE format which has to be indicated in 
 `configuration file <https://raw.githubusercontent.com/luidale/starpa/master/src/starpa/data/config.txt>`_ .
| File can not be sorted by position.

| BAM format is not currently supported.


- *identify*

| Aligned SE or pseudoSE reads in SAM format. 
| Reads require NH tag to describe the number of reported alignments.

| BAM format currently not supported.


- *cluster*

| Identified and quantified predicted processing products in BED format 
| (quantification at column #6).

|  folder bam:
| 	Aligned SE or pseudoSE reads in BAM format.
| 	Reads require NH tag to describe the number of reported alignments.

| If task "quantify" will be also executed:
| 	Additional input folder (given by parameter "quantify_sam_file_location"):
| 		Aligned SE or pseudoSE reads in SAM format 
| 		(BAM format currently not supported).
| 		Reads require NH tag to describe the number of reported alignments.


- *quantify*

| Predicted processing products in BED format (preferentially representatives form clustering).

| Additional input folder (given by parameter "quantify_sam_file_location"):
|	Aligned SE or pseudoSE reads in SAM format (BAM format currently not supported).
|	Reads require NH tag to describe the number of reported alignments.


**Output folder**

Output folder will contain parameter folder:

::

 parameters/
	eg. config.txt			-	copy of configuration file
	arguments.txt			-	command line arguments
	eg. libraries.txt		-	copy of library file
	eg. parameters.dev-2-100-2.txt	-	copy of Flaimapper-2 parameter file
 

Each task creates a subfolder with its name containing specific output 
of the task.

| XXX - library name
| strand - For or Rev
| Y -	order number of fragmented read group


- *trim*

::

 trim_info/
	XXX_triminfo.log	-	log of task
	XXX_triminfo.error	-	collected errors during trimming

 PE:
 discard/
	XXX_1_short.fq		-	forward reads discared while being too short after
					trimming
	XXX_2_short.fq		-	reverse reads discared while being too short after
					trimming
							
 XXX_trim_1.fq			-	trimmed forward reads
 XXX_trim_2.fq			-	trimmed reverse reads

 SE:
 discard/
	XXX_short.fq		-	reads discarded while being too short after 
					trimming
	XXX_untrimmed.fq	-	reads discarded while having no adapter trimmed
	
 XXX_trim.fq			-	trimmed reads

- *align*

::

 align_info/
	XXX_aligninfo.log	-	log of task
	
 XXX.sam			-	aligned reads

- *sam_sort*

::

 sort_info/
	XXX_sortinfo.log	-	log of task
	
 XXX_unmapped.sam		-	unmapped reads
 XXX_sort.sam			-	processed reads

- *pseudoSE*

::

 pseudoSE_info/
	XXX_pseudoSEinfo.log		-	log of task
	
 mismatched/
	XXX_pseudoSE_mismatch.sam	-	reads discarded while having too many
						mismatches
										
 too_many_matches/
	XXX_pseudoSE_multimatch.sam	-	reads discarded while haveing too many
						genomic matches
										
 XXX_pseudoSE.sam			-	processed reads
	
 If oligoA allowed:
 oligoA/
	XXX-oligoA-mm_pseudoSE.sam	-	reads with 3' oligoA (non-genome 
						encoded) which would have otherwise 
						discarded
	XXX-oligoA-pseudoSE.sam		-	reads with 3' oligoA (non-genome
						encoded)
	
- *identify*

::

 flaimapper/						
	flaimapper_info/
		XXX/
			XXX_strand_Y_flaimapper.information	-	log of flaimapper
			
	flaimapper_temp/
		XXX/
			XXX_strand_Y_flaimapper.tab		-	flaimapper predicitons
			
 bam/
	XXX_strand.bam						-	strand-wise sorted reads 
									from input
	XXX_strand.bam.bai					-	index of of bam file
	XXX_strand.sam 						-	NOT NEEDED
	
 identify_info/
	 XXX_strand_identifyinfo.log				-	log of task
	 
 XXX_strand_pp.BED						-	NOT NEEDED
 XXX_strand_pp_counted.BED					-	predicted processing 
									products with 
									quantification

			
- *cluster*

::

 cd_hit_est/
	pp_cd_hit_est.info		-	log of sequence identity based clustering 
						of combined and overlap clustered predicted
						processing products via CD-HIT-EST
	pp_combined.cdhit		-	genomic sequence of combined and overlap 
						clustered predicted processing products
	pp_combined.cdhit.clstr		-	clusters of combined and overlap clustered
						predicted processing products created via
						CD-HIT-EST
									
 contigs/
	XXX_contigs.BED			-	list of contigs identified
	XXX/
		contig_name.fasta	-	sequences of all reads belonging to the
						corresponding contigs
		contig_name.sam		-	all reads belonging to the
						corresponding contigs
									
 contigs_meta/
	combined_contigs_meta.BED	-	combined contigs to be used to create 
						metacontigs from all libraries
	XXX_contigs_meta.BED		-	list of contigs to be used to created
						metacontigs
	metacontig_cd_hit_est.info	-	log of sequence identity based clustering 
						of metacontigs via CD-HIT-EST
	metacontigs.cdhit		-	genomic sequence of metacontigs
	metacontigs.cdhit.clstr		-	clusters of metacontigs created via
						CD-HIT-EST
	metacontigs.BED			-	list of metacontigs in bed format
	pp_to_metacontig.BED		-	combined and overlap clustered predicted
						processing product match with metacontigs
						in BED-like format
									
 mpileup/
	XXX_strand_mpileup.info		-	log of bedtools mpileup
	
 wig/
	XXX_strand.wig			-	strand specific absolute read coverage
	XXX_strand_RPM.wig		-	strand specific relative read coverage
						as read per million mapped reads (RPM)
									
 pp_clusterinfo.log			-	log of task
 pp_unique.library_info			-	combined predicted processing 
						products and the origins of libraries
 pp_combined.BED			-	representatives of combined and overlap 
						clustered predicted processing products 
						in BED format
 pp_combined.cluster			-	overlap clusters of combined predicted 
						processing products
 pp_combined.library_info		-	representatives of combined and overlap 
						clustered predicted processing 
						products and the origins of libraries
 pp_metacontig.BED			-	representatives of predicted processing
						products from pp_combined.BED clustered
						by sequence identity supported by 
						metacontig clustering in BED format
 pp_metacontig.cluster			-	sequence identity clusters of predicted 
						processing products from pp_combined.BED
						supported by metacontig clustering

- *quantify*

::

 libraries/					-	data in library wise
	XXX.biotype_annotation.statistics	-	read alignement statistics
							by annotation biotypes
	XXX.gene_annotation.statistics		-	read alignement statistics
							by genes
	pp_metacontig_XXX_counted.BED		-	absolute quantification of 
							predicted processing products 
							in BED format
													
 collected.annotation2.statistics 		-	combined alignement	statistics
							by annotation biotypes
 pp_metacontig_biotype.BED			-	predicted processing products
							with biotype in BED-like format
 pp_metacontig_biotype_match.BED		-	predicted processing products
							match with genes in BED-like 
							format
 pp_metacontig_counts_total.tsv			-	absolute quantification of 
							predicted processing products 
							in BED format
 pp_metacontig_counts_RPM.tsv			-	relative quantification of 
							predicted processing products
							as read per million mapped reads
							(RPM) in BED format
 pp_metacontig_counts_biotype_RPM.tsv		-	relative quantification of 
							predicted processing products
							as RPM of biotype in BED format
 pp_metacontig_counts_groupped_biotype_RPM.tsv	-	relative quantification of 
							predicted processing products
							as RPM of biotype groups in BED 
							format
 pp_metacontig_cons_qual.tsv			-	quality of consensus sequence 
 							of predicted processing products
							expressed as frequency of the most
							abundant base in a given position
 pp_metacontig_cons_seq.tsv			-	consensus sequence of predicted 
							processing products
 pp_metacontig_coverage.tsv			-	coverage of reads assigned to 
							predicted processing products 
							at single position level
 pp_metacontig_genomic_seq.tsv			-	genomic sequence of predicted 
							processing products 
 pp_metacontig_rel_cov.tsv			-	relative coverage of predicted 
							processing products
 pp_metacontig_uniqness.tsv			-	mean number of genomic genomic 
							matches of reads assigned
							to the predicted processing 
							products

To do
-------------

Licence
-------
`GNU General Public License v3.0 <https://github.com/luidale/starpa/blob/master/LICENSE>`_

Authors
-------
`starpa` was written by `Hannes Luidalepp <luidale@gmail.com>`_