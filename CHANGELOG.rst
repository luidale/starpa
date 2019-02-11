=======
Changes
=======

v0.3.1 (pre-release)
-----------------

* PseudeSE - sequence quality conversion a la Edgar & Flyvbjerg, 2015
* Identify - pp quantification via featureCounts
* Cluster

  a)bam file is parsed using pybam (converted to Pyhton3).
  
  b)Parsing is conducted contig at the time which results low memory consumption (especially large genomes).
  
  c)added option to choose are contig data files (sam, fasta) created.
* Quantify

  a)statistic is collected over various tasks
  
  b)pp-s are selected by coverage

v0.3.0 (2018-04-24)
-----------------

* Added tests
* Identify - input sam is split fractionated by size by parsing
* Quantify - updated
* SE support checked

v0.2.1 (2018-03-23)
-----------------

* minor changes


v0.2 (2018-03-22)
-----------------

* commandline option added


v0.1 (2018-03-21)
-----------------

* initial release