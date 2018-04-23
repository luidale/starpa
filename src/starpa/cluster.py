#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Makes contigs.
Makes contigs for metacontigs.
Makes wig files.
Combines pp-s from different libraries.
Creates metacontigs and clusters pp-s

To do:
1)How to set min_cov_meta. It too high then some similar contigs are not
 clustered. What to look while setting:
 a)how many clustered pp-s are eventually combined together
2)Would parsing the sam file be faster than current mpileup,create contig sam and contig fasta.
Perhaps combining it all into single parsing of the sam file can be faster?
Coverage will be stored in contig wise (less memory), reads are kept in memory
as long the contig requirements are not fulfilled. Then reads are saved into fasta and sam file.
New reads are added until contig ends.
Before coverage pp-s have to be read in. When contig overlaping pp is finished,
then relative coverage of pp is checked and it pp is passed as qualified if rel cov is good enough.

BE CAREFUL!!
1)Coverage is read in for the whole genome.
 Big genomes (eukaryotic) can use too much memory.

'''

import os
import subprocess
import bisect
import shutil
import copy
import sys
import math
from itertools import groupby
import multiprocessing as mp

from pyfaidx import Fasta

class cluster():
    def __init__(self,settings,first_task):
        self.settings = settings
        self.make_folder(settings)
        self.cluster(settings,first_task)
        self.test_errors(settings)

    def make_folder(self,settings):
        '''
        Make output folder for task if needed
        '''
        if not os.path.exists(os.path.join(settings["--output"],"cluster")):
            os.makedirs(os.path.join(settings["--output"],"cluster"))
        if not os.path.exists(os.path.join(settings["--output"],"cluster","wig")):
            os.makedirs(os.path.join(settings["--output"],"cluster","wig"))
        if not os.path.exists(os.path.join(settings["--output"],"cluster","contigs")):
            os.makedirs(os.path.join(settings["--output"],"cluster","contigs"))
        if not os.path.exists(os.path.join(settings["--output"],"cluster","contigs_meta")):
            os.makedirs(os.path.join(settings["--output"],"cluster","contigs_meta"))
        if not os.path.exists(os.path.join(settings["--output"],"cluster","cd_hit_est")):
            os.makedirs(os.path.join(settings["--output"],"cluster","cd_hit_est"))
        if not os.path.exists(os.path.join(settings["--output"],"cluster","mpileup")):
            os.makedirs(os.path.join(settings["--output"],"cluster","mpileup")) 
            
    def cluster(self,settings,first_task):
        '''
        Clusters pp-s from different libraries
        '''
        pp_dict = {}
        strand_list = ["For", "Rev"]

        #MAKE INFO FILES
        f_represent = open(os.path.join(settings["--output"],"cluster","pp_combined.cluster"),"w")
        f_info = open(os.path.join(settings["--output"],"cluster","pp_clusterinfo.log"),"w")


        #PROCESS READS
        print("\tProcess reads")
        pool = mp.Pool(processes=settings["CPUs"])
        results = {library:pool.apply_async(self.process_reads, \
                          args = (settings,library,strand_list,first_task)) \
                          for library in sorted(settings["libraries"])}
        pp_dict = {library:p.get() for library,p in results.items()}
        pool.close()
        pool.join()
                           
        #COMBINE PP-s FROM DIFFERENT LIBRARIES
        print("\tCombine PP-s from different libraries")
        combined_pp_dict = {}
        for library in pp_dict:
            for pp in pp_dict[library]:
                combined_pp_dict.setdefault(pp,[]).append(library)
        #write pp library relations to file
        self.write_pp_library_info(sorted(list(combined_pp_dict.keys())),os.path.join(\
            settings["--output"],"cluster","pp_unique.library_info"),combined_pp_dict)
                
        #FINDS OVERLAPING BETWEEN PP-s
        pp_overlap = self.overlapping_pps(settings,combined_pp_dict)
                                
        #GET REPRESENTATIVE PP-s
        pp_list_representatives = self.get_pp_representatives(\
                                    pp_overlap,combined_pp_dict,f_represent)
        
        #WRITE PP-S TO THE FILE
##        self.write_pp_data(sorted(combined_pp_dict),os.path.join(\
##            settings["--output"],"cluster","pp_combined.BED"))
        self.write_pp_data(pp_list_representatives,os.path.join(\
            settings["--output"],"cluster","pp_combined.BED"),combined_pp_dict)
        self.write_pp_library_info(pp_list_representatives,os.path.join(\
            settings["--output"],"cluster","pp_combined.library_info"),combined_pp_dict)

        #COMBINE CONTIGS FOR METACONTIGS
        self.metacontigs(settings)

        #CLUSTER PP-s BY SEQUENCE AND METACONTIGS
        pp_representatives = self.cluster_pps(settings)
        

        #WRITE CLUSTERED PP-S TO THE FILE
        print("\tWrite clustered pp representatives to file")                      
        #print(pp_clust_represent_file)
        f_pp_clust = open(os.path.join(settings["--output"],\
                                       "cluster","pp_metacontig_unsorted.BED"),"w")
        for pp_rep in sorted(pp_representatives):
            strand = pp_rep.split("_")[1][0]
            chrom = "_".join(pp_rep.split("_")[1:-2])[1:]
            start = pp_rep.split("_")[-2].lstrip("0")
            end = str(int(start)+int(pp_rep.split("_")[-1]))
            f_pp_clust.write("\t".join([chrom,start,end,pp_rep,"1",strand])+"\n")
        f_pp_clust.close()

        #sort combined file
        print("\tSort clustered pp representatives") 
        sort_command = (
            settings["bedtools_call"], "sort",
            "-i",os.path.join(settings["--output"],"cluster","pp_metacontig_unsorted.BED"),
            ">",os.path.join(settings["--output"],"cluster","pp_metacontig.BED")
            )
        os.system(" ".join(sort_command))
        os.remove(os.path.join(settings["--output"],"cluster","pp_metacontig_unsorted.BED"))

        #remone unnecessary files
        os.remove(os.path.join(settings["--output"],"cluster","pp_combined.fa"))
        os.remove(os.path.join(settings["--output"],"cluster","contigs_meta",
                               "metacontigs.fasta"))

        #WRITE STATISTICS
        f_info.write("Uniq processing products:\t\t\t"+str(len(combined_pp_dict))+"\n")
        f_info.write("Processing products clustered by overlap:\t"+\
                     str(len(pp_list_representatives))+"\n")
        f_info.write("Reduction by clustering step:\t\t\t"+\
                     str(round((len(combined_pp_dict)-len(pp_list_representatives))/\
                               len(combined_pp_dict),3))+"\n")
        f_info.write("Processing products clustered by metacontigs:\t\t"+\
                     str(len(pp_representatives))+"\n")
        f_info.write("Reduction by clustering step:\t\t\t"+\
                     str(round((len(pp_list_representatives)-len(pp_representatives))/\
                               len(pp_list_representatives),3))+"\n")

        f_represent.close()
        f_info.close()
        
    def process_reads(self,settings,library,strand_list,first_task):
        '''
        Processing all read of a library into to a coverage dictionary.
        Cverage is used to create contigs, wig-files and select pp-s.
        '''
        #MAKE FOLDERS if needed
        if not os.path.exists(os.path.join(settings["--output"],"cluster","contigs",library)):
            os.makedirs(os.path.join(settings["--output"],"cluster","contigs",library))
        
        #CREATE DICTIONARY WITH COVERAGE
        coverage_dic = {}
        for strand_name in strand_list:
            print("\t\t"+ library,strand_name,"Read in coverage")
            if strand_name == "For":
                strand = "+"
            else:
                strand = "-"
            coverage_dic[strand] = self.mpileup_to_dic(settings,library,strand_name,first_task)
            
        #GET TOTAL READ NUMBER
        total_reads = self.get_total_reads(settings,library,first_task)
               
        ##CREATE CONTIGS
        self.get_contigs(settings,coverage_dic,"contigs",library)
        self.create_contig_data_files(settings,library,"contigs",coverage_dic[strand],\
                                      first_task)
        ##CREATE CONTIGS FOR METACONTIGS
        self.get_contigs(settings,coverage_dic,"contigs_meta",library)
        ##CREATE WIG FILES
        self.get_wig_file(settings,coverage_dic["+"],library,"For",total_reads)
        self.get_wig_file(settings,coverage_dic["-"],library,"Rev",total_reads)
        
        #READ IN PROCECCING PRODUCTS
        pp_dict_lib = self.read_in_pps(settings,library,coverage_dic,strand_list,first_task)

        return pp_dict_lib

    def mpileup_to_dic(self,settings,library,strand_name,first_task):
        '''
        Fills dictionary with positions in coverage in 1-based cordination system
        '''
        #gets coverage of the whole library in 1-based coordination system
        if first_task == "cluster":
            input_file = os.path.join(settings["--input"],"bam",library+"_"+strand_name+".bam")
        else:
            input_file = os.path.join(settings["--output"],"identify",\
                                      "bam",library+"_"+strand_name+".bam")

        #index bam if index missing
        if not os.path.isfile(input_file+".bai"):
            samtools_index_command = (
                settings["samtools_call"], "index",
                "-@", str(settings["samtools_threads"]),
                input_file
                )

            os.system(" ".join(samtools_index_command))

        #getting chromosome names
        idxstats_command = (
                        settings["samtools_call"], "idxstats",
                        input_file
                        )
        
        idxstats_output = subprocess.Popen(idxstats_command, stdout=subprocess.PIPE,\
                                    universal_newlines=True)
        chromosomes = []
        for line in idxstats_output.stdout:
            if not line.strip().split("\t")[0] == "*":
                chromosomes.append(line.strip().split("\t")[:2])

        #run pileup
        #counts coverage by 100 000 at once.
        pos_dic = {}
        for chrom in chromosomes:
            pos_intervals = list(range(0,int(chrom[1]),100000))+[int(chrom[1])]
            regions = [chrom[0]+":"+ str(x+1)+"-"+str(pos_intervals[i+1]) for i,x in enumerate(pos_intervals[:-1])]
            pos_dic[chrom[0]] = {}
            for region in regions:
                pileup_command = (
                                settings["samtools_call"], "mpileup",
##                                "-@", str(settings["samtools_threads"]),
                                "-d", "100000000", "-Bsf",
                                settings["genome"], input_file,
                                "--ff", "200", "-r", region
                                )
                proc = subprocess.Popen(pileup_command, stdout=subprocess.PIPE,stderr=open(os.path.join(\
                                            settings["--output"],"cluster","mpileup",\
                                            library+"_"+strand_name+"_mpileup.info"),"w"),\
                                            universal_newlines=True)
                
                for line in proc.stdout:
                    line = line.strip().split("\t")
                    #adds list:[read number,newly starting reads]
                    pos_dic[chrom[0]][int(line[1])]=[int(line[3]),int(line[4].count("^"))]
        return pos_dic
            
    def get_total_reads(self,settings,library,first_task):
        '''
        Gets total read number from pseudoSE info file
        '''
        if first_task in {"cluster", "identify"}:
            total_reads_file = os.path.normpath(os.path.join(settings["--input"],\
                                                settings["cluster"]["cluster_pseudoSE_location"],\
                                                library+"_pseudoSEinfo.log"))
        else:
            total_reads_file = os.path.join(settings["--output"],"pseudoSE","pseudoSE_info",\
                                  library+"_pseudoSEinfo.log")
        try:
            with open(total_reads_file) as f_in:
                for line in f_in:
                    if line.startswith("Passed reads:"):
                        total_reads = line.strip().split("\t")[1]
                        break
            return total_reads
        
        except FileNotFoundError:
            sys.exit("No file for total reads file ("+total_reads_file+") for library "+library)
        
    def get_contigs(self,settings,coverage_dic,filename_pref,library):
        '''
        Create contigs in BED format. Coordination system 0-based, end exclusive
        '''
        print("\t\t"+library,"Assamble",filename_pref)
        
        if filename_pref == "contigs":
            min_contig_length = settings["cluster"]["cluster_min_contig_length"]
            min_contig_cov = settings["cluster"]["cluster_min_contig_cov"]
            min_contig_reads = settings["cluster"]["cluster_min_contig_reads"]
            out_contigs = open(os.path.join(settings["--output"],"cluster","contigs",\
                                        library+"_"+filename_pref+".BED"),"w")
        else:
            min_contig_length = settings["cluster"]["cluster_min_contig_length_meta"]
            min_contig_cov = settings["cluster"]["cluster_min_contig_cov_meta"]
            min_contig_reads = settings["cluster"]["cluster_min_contig_reads_meta"]           
            out_contigs = open(os.path.join(settings["--output"],"cluster","contigs_meta",\
                                        library+"_"+filename_pref+".BED"),"w")

        #writing header only to contig file, not to contigs for metacontigs
        if filename_pref == "contigs":
            out_contigs.write("track name=\""+settings["libraries"][library][0]+\
                              "-"+settings["libraries"][library][1]+\
                              " ["+library+"]\""+\
                              " description=\""+settings["libraries"][library][0]+\
                              "-"+settings["libraries"][library][1]+\
                              " ["+library+"]\""+\
                              " useScore=1 visibility=pack color=100,50,0\n")
        for strand in sorted(coverage_dic):
            for chrom in sorted(coverage_dic[strand]):
                reads = 0
                #find the beginning of first suitable processing product 
                positions = list(sorted(coverage_dic[strand][chrom]))
                i = 0
                while i < len(positions):
                    if coverage_dic[strand][chrom][positions[i]][0] >= min_contig_cov:
                        start = positions[i]
                        prev_pos = positions[i]
                        reads = coverage_dic[strand][chrom][positions[i]][0]
                        i += 1
                        break                    
                    else:
                        i += 1
                            
                #continue processing positions
                for pos in positions[i:]: 
                    #contig starts when new chrom or position bigger than previous+1
                    if pos != prev_pos+1:
                        #testing minimum read number of the contig
                        if reads >= min_contig_reads:
                            #testing minimum contig length
                            if (prev_pos-start+1) >= min_contig_length:
                                #write contig to file
                                name = chrom+strand+str(start-1)+"_"+str(prev_pos)+"_"+str(reads)
                                out_contigs.write("\t".join(\
                                    [chrom,str(start-1),str(prev_pos),name,"1",strand+"\n"]))
                        #start new contig if minimum coverage at position is fulfilled
                        if coverage_dic[strand][chrom][pos][0] >= min_contig_cov:
                            prev_pos = pos
                            start = pos
                            reads = coverage_dic[strand][chrom][pos][0] #line[4].count("^")
                        else:
                            reads = 0
                                      
                    else:
                        #continueing of existing contig
                        #if minimum coverage at position is fulfilled
                        if coverage_dic[strand][chrom][pos][0] >= min_contig_cov:
                            prev_pos = pos
                            #adds number of reads starting at this position
                            reads += coverage_dic[strand][chrom][pos][1] 

                #write last contig to the file
                #testing minimum read number of the contig
                if reads >= min_contig_reads:
                    #testing minimum contig length
                    if (pos-start+1) >= min_contig_length:
                        #write contig to file
                        name = chrom+strand+str(start-1)+"_"+str(prev_pos)+"_"+str(reads)
                        out_contigs.write("\t".join([chrom,str(start-1),\
                                                     str(prev_pos),name,"1",strand+"\n"]))
        out_contigs.close()

    def create_contig_data_files(self,settings,library,filename_pref,coverage_dic,\
                                 first_task):
        '''
        Make sam and fasta files for the contigs.
        '''
        #take contigs from BED file and creates sam, fasta
        contig_output_folder = os.path.join(settings["--output"],"cluster","contigs",library)
        if not os.path.exists(contig_output_folder):
            os.makedirs(contig_output_folder)

        #input folder
        if first_task == "cluster":
            input_folder = settings["--input"]
        else:
            input_folder = os.path.join(settings["--output"],"identify")

        #index bam if index missing
        for strand_name in ["For","Rev"]:
            bam_file = os.path.join(input_folder,"bam",\
                                 library+"_"+strand_name+".bam")
                
            if not os.path.isfile(bam_file+".bai"):
                samtools_index_command = (
                    settings["samtools_call"], "index",
                    "-@", str(settings["samtools_threads"]),
                    bam_file
                    )

                os.system(" ".join(samtools_index_command))

        #read in contigs    
        with open(os.path.join(settings["--output"],"cluster","contigs",\
                               library+"_"+filename_pref+".BED")) as f_in:

            #remove header
            f_in.readline()
                        
            for line in f_in:
                line = line.strip().split("\t")
                chrom,start,end,strand,name = line[0],line[1],line[2],line[5],line[3]
                if strand == "-":
                    strand_name = "Rev"
                else:
                    strand_name = "For"
                    
                #create SAM file
                contig_reads_command = (
                        settings["samtools_call"], "view",
                        "-@", str(settings["samtools_threads"]-1),
                        os.path.join(input_folder,"bam",\
                                     library+"_"+strand_name+".bam"),
                        '"'+chrom+":"+str(start)+"-"+str(end)+'"',
                        "-o", os.path.join(contig_output_folder,name+".sam")
                        )
                os.system(" ".join(contig_reads_command))

                #create FASTA file
                contig_fasta_command = (
                    settings["samtools_call"], "view",
                    "-@", str(settings["samtools_threads"]-1),
                    os.path.join(input_folder,"bam",\
                                     library+"_"+strand_name+".bam"),
                    chrom+":"+str(start)+"-"+str(end),"|",
                    "awk", '\'{OFS="\\t"; print ">"$1"\\n"$10}\'',
                    ">", os.path.join(contig_output_folder,name+".fasta")
                    )
                os.system(" ".join(contig_fasta_command))
                        
                ##if contig is in reverse strand,
                ##the fasta sequence have to be reverse complemented
                if strand == "-":
                    #create reverse complement file
                    f_out = open(os.path.join(contig_output_folder,name+".fasta_temp"),"w")
                    with open(os.path.join(contig_output_folder,name+".fasta")) as f_in:
                        for line in f_in:
                            f_out.write(line)
                            f_out.write(self.rev_comp4(f_in.readline().strip())+"\n")
                    f_out.close()
                    #delete original file and rename temporary file
                    os.remove(os.path.join(contig_output_folder,name+".fasta"))
                    os.rename(os.path.join(contig_output_folder,name+".fasta_temp"),\
                              os.path.join(contig_output_folder,name+".fasta"))

##                #create HTML file
##                self.create_contig_html(settings,chrom,start,end,strand,name,coverage_dic,total_reads)       

    def rev_comp4(self,seq):
        '''
        Converts DNA sequence to reverse complement
        '''
        # too lazy to construct the dictionary manually, use a dict comprehension
        seq1 = 'ATCGNRYMKBDHV *TAGCNYRKMVHDB *'
        seq_dict = { seq1[i]:seq1[i+15] for i in range(30) if i < 15}
        return "".join([seq_dict[base] for base in reversed(seq)])
                
    def create_contig_html(self,settings,chrom,start,end,strand,name,coverage_dic,total_reads):
        '''
        Creates html file for contigs
        '''
        #FUNCTION IS NOT READY
        #COLLECTION OF DATA FOR HTML IS MISSING
        #This should be done while parsing pileup
        #Perhaps should wait until pp-s are annotated and then add those alo here.
        with open(os.path.join(settings["--output"],"cluster",library,name+".sam")) as f_sam:
            nr_reads = 0
            genomic_uniqness = 0
            for line in f_sam:
                line=line.strip().split("\t")
                nr_reads += 1
                genomic_uniqness += int(line[-1].strip("NH:i:"))
            genomic_unicness = round(genomic_uniqness/nr_reads,2)
            nr_reads_rpm = round(nr_reads/total_reads,3)
            pp_positions = list(range(start,end)) 
            max_cov = max([int(coverage_dic[strand][chrom][x+1][0]) for x in pp_positions])
    
            

        #with open(os.path.join(settings["--output"],"cluster",library,name+".html"),"w") as f_html:
            #f_html.write("<HTML>\n<H2>"+name"</H2>\n")
            #f_html.write("<P>Length: "+str(int(end)-int(start)+"\n")
            #f_html.write("<P>Number of reads aligned: "+str(nr_reads)+"\n")
            #f_html.write("<P>Normalized expression level (RPM): "+str(nr_reads_rpm)+"\n")
            #f_html.write("<P>Maximum coverage: "+str(max_cov)+"\n")
            #f_html.write("<P>Genomic uniqness: "+str(genomic_unicness)+"\n")
            #f_html.write("<P>Genomic position: "+chrom+":"+str(start)+"-"+str(end)+"\n")
##            #f_html.write("<P>Genome mapped repeats overlapping the loci: ????\n")
##            #f_html.write("<P>Alignment of read-derived sequence consensus to genomic sequence:\n<P>")
##            #f_html.write("<font face=\"\'Courier New\', monospace\">Genomic:&nbsp;&nbsp;&nbsp;$"+genomic_seq+"<BR>\n")
##            #f_html.write("Consensus:&nbsp;"+$cons_seq+"<BR>\n")
##            #f_html.write("Cons.qual&nbsp;&nbsp;"+$seq_quality{+"</font>\n")
##            #f_html.write(
##            #f_html.write("<P>Annotation features overlapping the contig:\n<P><TABLE border=1>
##            #f_html.write("<TR><td><B>Feature ID</B></td><td><B>Exon/Intron</B></td><td><B>Orientation</B></td><td><B>Feature status</B></td><td><B>Feature biotype</B></td><td><B>Feature alternative name</B></td><td><B>Feature description</B></td></TR>\n";
##            #f_html.write(
##            #f_html.write("<P>Other genomic loci with partial sequence match: NA\n")
            #f_html.write("<P><A href=\""+os.path.join(settings["--output"],"cluster",library,name+".fasta")+"\">Download sequences of the reads</A>\n<P>")
            #f_html.write("<P><A href=\""+os.path.join(settings["--output"],"cluster",library,name+".sam")+"\">Download original alignment of the reads in SAM format</A>\n<P>")

    def get_wig_file(self,settings,coverage_dic,library,strand_name,total_reads):
        '''
        Creates wig files
        '''
        print("\t\t"+library,strand_name,"Create wig files")
        wig = open(os.path.join(settings["--output"],"cluster","wig",\
                                library+"_"+strand_name+".wig"),"w")
        wig_rpm = open(os.path.join(settings["--output"],"cluster","wig",\
                                    library+"_"+strand_name+"_RPM.wig"),"w")
        ##ADD TRACK DEFINITION LINE
        track_line =["track type=wiggle_0"]
        track_line_rpm =["track type=wiggle_0"]
        #add name
        track_line.append("name=\""+settings["libraries"][library][0]+\
                          "-"+settings["libraries"][library][1]+\
                          " "+strand_name+" ["+library+"]\"")
        track_line_rpm.append("name=\""+settings["libraries"][library][0]+\
                              "-"+settings["libraries"][library][1]+\
                              " "+strand_name+" ["+library+"] RPM\"")
        #description
        track_line.append("description=\""+settings["libraries"][library][0]+\
                          "-"+settings["libraries"][library][1]+\
                          " "+strand_name+" ["+library+"]\"")
        track_line +=["visibility=full","autoScale=on","alwaysZero=on",\
                      "maxHeightPixels=50:50:5"]
        track_line_rpm.append("description=\""+settings["libraries"][library][0]+\
                              "-"+settings["libraries"][library][1]+\
                              " "+strand_name+" ["+library+"] RPM\"")
        track_line_rpm +=["visibility=full","autoScale=on","alwaysZero=on",\
                          "maxHeightPixels=50:50:5"]
        #add color
        if strand_name == "For":
            track_line.append("color=0,120,0\n")
            track_line_rpm.append("color=0,120,0\n")
            orientation =""
        else:
            track_line.append("color=0,0,255\n")
            track_line_rpm.append("color=0,0,255\n")
            orientation ="-"
        #print(track_line)
        wig.write(" ".join(track_line))
        wig_rpm.write(" ".join(track_line_rpm))

        #ADD COVERAGE DATA
        for chrom in sorted(coverage_dic):
            prev_pos = -1
            for pos in sorted(coverage_dic[chrom]):
                if pos > prev_pos +1:
                    wig.write("variableStep chrom="+chrom+"\n")
                    wig_rpm.write("variableStep chrom="+chrom+"\n")
                wig.write("\t".join([str(pos),orientation+str(coverage_dic[chrom][pos][0])+"\n"]))
                wig_rpm.write("\t".join([str(pos),orientation+\
                                         str(round(coverage_dic[chrom][pos][0]/\
                                                   float(total_reads)*1000000,5))+"\n"]))
                prev_pos = pos
        wig.close()
        wig_rpm.close()

    def read_in_pps(self,settings,library,coverage_dic,strand_list,first_task):
        '''
        Reads in pp-s from file and selects them according to the read number
        and relative coverage.
        '''
        print("\t\t"+library,"Read in processing products")
        pp_list = []
        for strand_name in strand_list:
            if first_task == "cluster":
                input_BED = open(os.path.join(settings["--input"],\
                                            library+"_"+strand_name+\
                                            settings["cluster"]["cluster_input_file_suffix"]))
            else:
                input_BED = open(os.path.join(settings["--output"],"identify",\
                                            library+"_"+strand_name+"_pp_counted.BED"))
            for line in input_BED:
                while line[0] == "#": #remove header
                    line = input_BED.readline()
                line = line.strip().split("\t")
                #CONSIDER PP-S ONLY WITH ENOUGH READS
                if int(line[6]) < settings["min_pp_reads"]:
                    continue
                chrom = line[0]
                start = int(line[1])
                end = int(line[2])
                name = line[3]
                strand = line[5]
                count = int(line[6])
                #CONDSIDER PP-S WITH REASONABLE RELATIVE COVERAGE
                #pp_coverage is average over pp.
                #Could use also otherthings (eg. coverage at ends).
                pp_positions = list(range(start,end)) 
                pp_coverage = sum([int(coverage_dic[strand][chrom][x+1][0]) \
                                   for x in pp_positions])/len(pp_positions)
                ##gets index for relative coverage
                ##in short: it will compares length of the pp with the size range
                ##and gives index for that
                index = bisect.bisect(settings["cluster"]["cluster_rel_cov_size_range"], len(pp_positions))
                if count < pp_coverage*settings["cluster"]["cluster_rel_cov_list"][index]:
                    continue
                pp_list.append((chrom,start,end,strand,name))
            input_BED.close()
        return pp_list

    def overlapping_pps(self,settings,combined_pp_dict):
        '''
        Creates dictionary with pp overlapping.
        Value for each pp is list of pp-s it is overlapping
        '''
        print("\tFinds overlapping processing products")
        pp_overlap = {} #dictionary of pp-s with overlapping pp-s
        #goes through all processing products one by one
        for i, pp in enumerate(sorted(combined_pp_dict)):
            #adds pp to the dictionary
            if pp not in pp_overlap:
                pp_overlap[pp] = []
            #tests the processing product against all following
            #[as preceedings are already tested] processing products
            for pp_test in sorted(combined_pp_dict)[i+1:]:
                #skip if they are not in same chromosome
                if pp[0] != pp_test[0]:
                    break
                #skip if they are not in same strand
                if pp[3] != pp_test[3]:
                    continue
                #skip if the pp_test is too far in upstream it is excluded
                if pp[1]-pp_test[1] > settings["max_length"]*2:
                    continue
                #if the pp_test is too far in downstream the testing for the pp is stopped
                if pp_test[2]-pp[2] > settings["max_length"]*2:
                    break        
                #testing the overlap
                ##gives length of overlap, BED coordination system checked
                overlap_len = len(range(max(pp[1],pp_test[1]),min(pp[2],pp_test[2]))) 
                if overlap_len != 0:
                    #testing overlap length for both (initial and tested) pps
                    ##enough overlap from length of initial pp
                    if (pp[2]-pp[1]-settings["non_overlap"]) <= overlap_len:
                        ##enough overlap from length of tested
                        if (pp_test[2]-pp_test[1]-settings["non_overlap"]) <= overlap_len:
                            pp_overlap[pp].append(pp_test)
                            #adds tested pp to the overlap list of the pp
                            if pp_test not in pp_overlap: 
                                pp_overlap[pp_test]=[pp]
                            else: #adds inital to overlap list of the tested pp
                                pp_overlap[pp_test].append(pp)
        return pp_overlap

    def get_pp_representatives(self,pp_overlap,combined_pp_dict,f_represent):
        '''
        Creates list of representative pp-s and writes them and theim with the
        represented pp-s to the file.
        '''
        print("\tGet representative processing products")
        checked_pp = set() # checked pp-s
        pp_list_representatives = [] #selected pp-s which are unique.
        for pp in sorted(pp_overlap.keys()):
            ##Get groups of overlapping processing products
            #if pp is checked earlier it will not be checked again (it will be skipped)
            if pp in checked_pp:
                continue
            #pp-s which overlap with something
            if len(pp_overlap[pp])!= 0:
                #getting pp group where overlaping pp-s are
                #connected with some overlaping pp-s with each other.
                pp_group = self.get_pp_group({},pp,pp_overlap)
                checked_pp.update(set(pp_group))
                
                ##Getting representativ pp-s for each group
                #Strategy:
                #Take pp which overlaps with the highest number
                #(if several with equal coverage take with medium middle genomic position).
                #Then recalculate the overlap of the rest by removing representative
                #and pp it is overlaping and take again highest.
                #Repeat until all are covered
                #group_representative_pp = self.get_representatives_overlaps(pp_group,f_represent,pp[3])
                group_representative_pp = self.get_representatives_libraries(pp_group,\
                                                                             combined_pp_dict,\
                                                                             f_represent,pp[3])
                pp_list_representatives += group_representative_pp
            #pp-s which do not overlap with anything are added to the reduced list
            elif len(pp_overlap[pp])== 0:
                f_represent.write(pp[4]+"\n")
                pp_list_representatives.append(pp)
        pp_list_representatives = sorted(pp_list_representatives)
        return pp_list_representatives

    def metacontigs(self,settings):
        '''
        Creates metacontigs and clusters them
        '''
        with open(os.path.join(settings["--output"],"cluster","contigs_meta",\
                               "combined_contigs_meta_unsorted.BED"),'wb') as wfd:
            for library in settings["libraries"]:
                with open(os.path.join(settings["--output"],"cluster","contigs_meta",\
                                       library+"_contigs_meta.BED"),'rb') as fd:
                    shutil.copyfileobj(fd, wfd, 1024*1024*10)

        with open(os.path.join(settings["--output"],"cluster","contigs_meta",\
                               "combined_contigs_meta_unsorted.BED")) as f_in:
            for line in f_in:
                print("a",line)
                
        #sort combined file
        sort_command = (
            settings["bedtools_call"], "sort",
            "-i",os.path.join(settings["--output"],"cluster","contigs_meta",\
                              "combined_contigs_meta_unsorted.BED"),
            ">",os.path.join(settings["--output"],"cluster","contigs_meta",\
                             "combined_contigs_meta.BED")
            )
        os.system(" ".join(sort_command))
        os.remove(os.path.join(settings["--output"],"cluster","contigs_meta",\
                               "combined_contigs_meta_unsorted.BED"))

        with open(os.path.join(settings["--output"],"cluster","contigs_meta",\
                              "combined_contigs_meta.BED")) as f_in:
            for line in f_in:
                print("b",line)
        
        #merge contigs into metacontigs
        merge_command = (
            settings["bedtools_call"], "merge",
            "-i",os.path.join(settings["--output"],"cluster","contigs_meta",\
                              "combined_contigs_meta.BED"),
            "-s", "-c", "5,6", "-o", "disttinct,distinct",
            ">", os.path.join(settings["--output"],"cluster","contigs_meta",\
                              "metacontigs.BED")
            )
        os.system(" ".join(merge_command))
        
        with open(os.path.join(settings["--output"],"cluster","contigs_meta",\
                              "metacontigs.BED")) as f_in:
            for line in f_in:
                print("c",line)
                
        #name metacontigs
        self.name_contigs(settings,os.path.join(settings["--output"],"cluster","contigs_meta",\
                                                "metacontigs.BED"))

        with open(os.path.join(settings["--output"],"cluster","contigs_meta",\
                              "metacontigs.BED")) as f_in:
            for line in f_in:
                print("f",line)
                
        #create fasta for metacontigs
        fasta_command = (
            settings["bedtools_call"], "getfasta",
            "-fi", settings["genome"], "-s",
            "-bed", os.path.join(settings["--output"],"cluster","contigs_meta",\
                                 "metacontigs.BED"),
            ">", os.path.join(settings["--output"],"cluster","contigs_meta",\
                              "metacontigs.fasta")
            )
        os.system(" ".join(fasta_command))

        with open(os.path.join(settings["--output"],"cluster","contigs_meta",\
                              "metacontigs.fasta")) as f_in:
            for line in f_in:
                print("d",line)


        self.rename_fasta(os.path.join(settings["--output"],"cluster","contigs_meta",\
                                       "metacontigs.fasta"),\
                     os.path.join(settings["--output"],"cluster","contigs_meta",\
                                  "metacontigs.BED")) 

        #cluster metacontigs
        print("\tCluster metacontigs")
        cluster_command = (
            settings["cd_hit_est_call"],
            "-i", os.path.join(settings["--output"],"cluster","contigs_meta",\
                               "metacontigs.fasta"),
            "-o", os.path.join(settings["--output"],"cluster","contigs_meta",\
                               "metacontigs.cdhit"),
            "-M", "3000", "-d", "0", "-r", "0", "-c", "0.9", "-s", "0.5",
            ">", os.path.join(settings["--output"],"cluster","contigs_meta","metacontig_cd_hit_est.info")
            )
            #-c - sequence identity threshold
            #-s - length difference cutoff
        os.system(" ".join(cluster_command))
        
    def cluster_pps(self,settings):
        '''
        Clusters pp-s according to the metacontigs
        '''
        #Write pp-BEDfile to fasta file
        print("\tWrite pp-file to fasta format")
        fasta_command = (
            settings["bedtools_call"], "getfasta",
            "-fi", settings["genome"], "-s", 
            "-bed", os.path.join(settings["--output"],"cluster","pp_combined.BED"),
            ">", os.path.join(settings["--output"],"cluster","pp_combined.fa")
            )
        os.system(" ".join(fasta_command))

        self.rename_fasta(os.path.join(settings["--output"],"cluster","pp_combined.fa"),\
                          os.path.join(settings["--output"],"cluster","pp_combined.BED"))

        #Cluster pp-s with CD-HIT-EST
        print("\tCluster processing products by sequence")
        cluster_pp_command =(
            settings["cd_hit_est_call"],
            "-i", os.path.join(settings["--output"],"cluster","pp_combined.fa"),
            "-o", os.path.join(settings["--output"],"cluster","cd_hit_est",\
                               "pp_combined.cdhit"),
            "-M", "3000","-d", "0", "-r", "0", "-c", "1","-s","1",
            ">", os.path.join(settings["--output"],"cluster","cd_hit_est","pp_cd_hit_est.info")
            )
        os.system(" ".join(cluster_pp_command))

        #Generate dictionary with metacontig relations
        print("\tGet metacontig relations")
        metacontig_dict_relation = self.get_feature_rel_dic(\
            os.path.join(settings["--output"],"cluster","contigs_meta",\
                         "metacontigs.cdhit.clstr"),True)

        #Generate dictionary with pp relations
        print("\tGet pp relations")
        pp_dict_relation = self.get_feature_rel_dic(\
            os.path.join(settings["--output"],"cluster","cd_hit_est",\
                         "pp_combined.cdhit.clstr"),False)
                
        ##Check are clustered pp-s in same metacontig group

        #Find for each pp a metacontig match
        print("\tFind metacontig matches for pp-s")
        intersect_command = (
            settings["bedtools_call"], "intersect",
            "-a", os.path.join(settings["--output"],"cluster","pp_combined.BED"),
            "-b", os.path.join(settings["--output"],"cluster","contigs_meta",\
                               "metacontigs.BED"),
            "-wao","-s",
            ">", os.path.join(settings["--output"],"cluster","contigs_meta",\
                              "pp_to_metacontig.BED")
            )
        os.system(" ".join(intersect_command))

        #Read in pp-metacontig matches
        f_pp_metacont = open(os.path.join(settings["--output"],"cluster","contigs_meta",\
                                          "pp_to_metacontig.BED"))
        pp_metacont_dic = {}
        for line in f_pp_metacont:
            pp = line.strip().split("\t")[3]
            metacont = line.strip().split("\t")[9]
            #Wrong?:#list as there could be two metacontigs overlaping the pp
            pp_metacont_dic[pp]=metacont 
        f_pp_metacont.close()
        
        #Verify pp sequence based clustering with metacontig clustering
        pp_representatives = {} #list of representative pp-s from clustering by sequence
        pp_groupped = set() #list of pp which are represented by representative
        for pp in sorted(pp_dict_relation):
                #goes through all pp-s in the group and looks up metacontigs
                #first pp-s with uniq metacontigs are considered as representatives
                #eg. first pp in group will be all cases representative
                #also next pp-s which have does not have metacontig
                #which is the same as representatives of the group
            pp_representatives[pp] = []
            group_representative = [pp]

            for pp_in_group in pp_dict_relation[pp][1:]:
                #test does metacontig of pp overlap with metacontigs of group representatives
                new_group_representative = copy.deepcopy(group_representative)
                for representative in new_group_representative:
                    #if metacontig clustering then it is groupped
                    if pp_metacont_dic[representative] in \
                            metacontig_dict_relation[pp_metacont_dic[pp_in_group]]:
                        pp_representatives[representative].append(pp_in_group)
                        break
                #no metacontig clustering with any group representatives
                #new representative
                else:
                    pp_representatives[pp_in_group] = []
                    group_representative.append(pp_in_group)



        #Write clustering file
        f_out = open(os.path.join(settings["--output"],"cluster","pp_metacontig.cluster"),"w")
        for pp in sorted(pp_representatives):
            if pp_representatives[pp] == []:
                f_out.write(pp+"\n")
            else:
                f_out.write("\t".join([pp,"-"]+\
                                    [";".join([x for x in pp_representatives[pp]])])+"\n")              
        return pp_representatives     
    
    def get_pp_group(self,pp_group,pp,pp_overlap):
        '''
        Recursive function which will collects a pp group where pp-s are overlaping
        '''
        #adds new subgroup
        if pp not in pp_group:
            pp_group[pp] = pp_overlap[pp]
            for pp_new in pp_group[pp]: #collects all overlaps of overlaping pp- of initial pp
                pp_group = self.get_pp_group(pp_group,pp_new,pp_overlap)
        return pp_group

    def get_representatives_overlaps(self,pp_group,f_represent,strand):
        '''
        Function to get representative pp-s.
        Pp having highest number of overlaps.
        '''
        representatives = []
        recalc_pp_group = copy.deepcopy(pp_group)
        while len(recalc_pp_group) != 0:
            max_overlap_pp = self.get_max_overlaps_pp(recalc_pp_group,f_represent,strand)
            representatives.append(max_overlap_pp)
            recalc_pp_group = self.recalc_pp_overlap(recalc_pp_group,max_overlap_pp)
        return representatives

    def get_representatives_libraries(self,pp_group,combined_pp_dict,f_represent,strand):
        '''
        Function to get representative pp-s.
        Pp found in highest number of libraries.
        '''
        representatives = []
        recalc_pp_group = copy.deepcopy(pp_group)
        while len(recalc_pp_group) != 0:
            max_overlap_pp = self.get_max_libraries_pp(recalc_pp_group,combined_pp_dict,\
                                                       f_represent,strand)
            representatives.append(max_overlap_pp)
            recalc_pp_group = self.recalc_pp_overlap(recalc_pp_group,max_overlap_pp)
        return representatives

    def get_max_overlaps_pp(self,pp_group,f_represent,strand):
        '''
        Finds processing product with maximum number of overlaps
        '''
        #highest number of overlapping pp-s:
        highest=max([len(pp_group[x]) for x in pp_group])
        #creates list of pp-s with highest overlapping pp-s:
        max_overlap_pp_list = [k for k,v in pp_group.items() if len(v) == highest] 
        list_len=len(max_overlap_pp_list)
        #sorts pp-s by location of middle position
        max_overlap_pp_list=sorted(max_overlap_pp_list, key= lambda x: x[1] +x[2])
        if strand == "+":
            representative = max_overlap_pp_list[math.floor(list_len/2)]
        else:
            representative = max_overlap_pp_list[math.ceil(list_len/2)-1]
        f_represent.write("\t".join([representative[4],"-"]+\
                                    [";".join([x[4] for x in pp_group[representative]])])+"\n")
        return representative #takes middle pp in the list

    def get_max_libraries_pp(self,pp_group,combined_pp_dict,f_represent,strand):
        '''
        Finds processing product with maximum number of libraries
        '''
        #highest number of libraries:
        highest=max([len(combined_pp_dict[x]) for x in pp_group])
        #creates list of pp-s with highest libraries:
        max_overlap_pp_list = [k for k,v in pp_group.items() if len(combined_pp_dict[k]) == highest]
        list_len=len(max_overlap_pp_list)
        #sorts pp-s by location of middle position:
        max_overlap_pp_list=sorted(max_overlap_pp_list, key= lambda x: x[1] +x[2]) 
        if strand == "+":
            representative = max_overlap_pp_list[math.floor(list_len/2)]
        else:
            representative = max_overlap_pp_list[math.ceil(list_len/2)-1]
        f_represent.write("\t".join([representative[4],"-"]+\
                                    [";".join([x[4] for x in pp_group[representative]])])+"\n")
        return representative #takes middle pp in the list

    def recalc_pp_overlap(self,pp_group,max_overlap_pp):
        '''
        Removes represented pp-s from the list.
        '''
        reduced_pp_group = {}
        for pp in pp_group:
            #removes max_overlap_pp and its overlaps
            if pp in pp_group[max_overlap_pp] or pp==max_overlap_pp:
                continue
            #removes from overlap list max_overlap_pp and pp-s overlaping with max_overlap_pp
            reduced_overlap = list(set(pp_group[pp])-\
                                   set(pp_group[max_overlap_pp]).union({max_overlap_pp}))
            reduced_pp_group[pp] = reduced_overlap
        return reduced_pp_group
                           
    def write_pp_data(self,pp_list,file_name,combined_pp_dict):
        '''
        Writes pp-s in to BED file.
        '''
        #bed file
        f_out_BED = open(file_name,"w")
        f_out_BED.write("#Chr\tStart\tEnd\tGeneID\tScore\tStrand\n")
        for pp in pp_list:
            f_out_BED.write(pp[0]+"\t")
            f_out_BED.write(str(pp[1])+"\t")
            f_out_BED.write(str(pp[2])+"\t")
            f_out_BED.write(str(pp[4])+"\t")
            f_out_BED.write(str(0)+"\t")
            f_out_BED.write(pp[3]+"\n")
        f_out_BED.close()
        
    def write_pp_library_info(self,pp_list,file_name,combined_pp_dict):
        '''
        Write pp and library relations to file
        '''
        #pp-library origin file
        f_out_BED = open(file_name,"w")
        for pp in pp_list:
            f_out_BED.write("\t".join([pp[4],"-"]+\
                                      [";".join([library for library in combined_pp_dict[pp]])])+"\n")
        f_out_BED.close()
        
    def name_contigs(self,settings,input_file):
        '''
        Names contigs
        '''
        #get genome
        genome = Fasta(settings["genome"],one_based_attributes=False)
        #open input and output
        f_in = open(input_file)
        f_out = open(input_file+"2","w") #file with new name
        #parsing input and writing fasta of metacontig
        print("\tName metacontigs")
        for line in f_in:
            print("e",line)
            chrom = line.strip().split("\t")[0]
            start = int(line.strip().split("\t")[1])
            end = int(line.strip().split("\t")[2])
            strand = line.strip().split("\t")[3]
            score = line.strip().split("\t")[4]
            name = "MC_"+strand+chrom+"_"+\
                    (len(str(len(genome[chrom])))-len(str(start)))*"0"+\
                    str(start)+"_"+str(end-start+1)
            f_out.write("\t".join([chrom,str(start),str(end),name,score,strand])+"\n")
        f_in.close()
        f_out.close()
        #old file will be replaced by new file
        os.remove(input_file)
        os.rename(input_file+"2",input_file)

##    def genome_parser(self,fasta_name):
##        '''
##        given a fasta file. yield tuples of header, sequence
##        Genome seq starts from position 0
##        '''
##        fh = open(fasta_name)
##        # ditch the boolean (x[0]) and just keep the header or sequence since
##        # we know they alternate.
##        faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
##        for header in faiter:
##            # drop the ">"
##            headerStr = header.__next__()[1:].strip()
##            # join all sequence lines to one.
##            seq = "".join(s.strip() for s in faiter.__next__())
##            yield (headerStr, seq)
        
    def rename_fasta(self,input_file,name_file):
        '''
        Renames sequencies in fasta file
        '''
        f_in = open(input_file)
        f_out = open(input_file+"2","w")
        f_name = open(name_file)
        for line in f_name:
            while line[0] == "#":
                line = f_name.readline()
            name = line.strip().split("\t")[3]
            f_out.write(">"+name+"\n")
            f_in.readline()
            f_out.write(f_in.readline())
        f_in.close()
        f_out.close()
        f_name.close()
        os.remove(input_file) #remove input file
        os.rename(input_file+"2",input_file) #change name of new file

    def get_feature_rel_dic(self,feature_file,add_whole_group):
        '''
        Creates dictionary relating features
        '''
        feature_dic = {}
        f_feature = open(feature_file)
        f_feature.readline() #remove first line
        group = []
        for line in f_feature:
            #if not cluster header then feature is added to the group
            if line.find(">Cluster") == -1:
                feature = line.strip().split(">")[1].split("...")[0]
                group.append(feature)
            else:
                #if new group header, then adds each member of the group to the
                #dictionary as keys. other members are as list the values
                if add_whole_group: #all group members will be added as keys
                    for i,feature in enumerate(group):
                        feature_dic[feature] = group[:i]+group[i+1:]
                else: #only representative (first in the group) of the group will\
                      #be added as key and the whole group is as values
                    feature_dic[group[0]] = group
                group=[] #clean group
                
        #add last feature group
        if add_whole_group: #all group members will be added as keys
            for i,feature in enumerate(group):
                feature_dic[feature] = group[:i]+group[i+1:]
        else: #only representative (first in the group) of the group will\
              #be added as key and the whole group is as values
            feature_dic[group[0]] = group
        f_feature.close()
        return feature_dic

    def test_errors(self,settings):
        '''
        Check errors
        '''
        info_file = os.path.join(settings["--output"],"cluster","pp_clusterinfo.log")
        if not os.path.isfile(info_file):
            sys.exit('Task "cluster" incomplete. Infofile ' + info_file + " missing")
        wanted_read_beginnings = {"Uniq","Processing","Reduction"}
        with open(info_file) as f_in:
            read_beginnings = set()
            for line in f_in:
                read_beginnings.add(line.split(" ")[0].strip())
            if not wanted_read_beginnings.issubset(read_beginnings):
                sys.exit('Error in task "cluster", check infofile: '+info_file)
