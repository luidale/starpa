#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Converts sam file to bam file
Identifies processing products (pp) via flaimapper
Converts Flaimapper output to BED format.
Counts reads for pp-s

TO DO:
1)Is which counting is fastest:
 a)Current: intersect via fragmented bed
 b)Parse through whole sam file (as in quantify=
 c)Featurecounts when it will have parameter to set non-overlapped bases

'''

from itertools import groupby
import os
import bisect
import shutil
import sys
import multiprocessing as mp

from pyfaidx import Fasta

class identify():
    def __init__(self,settings,first_task):
        self.settings = settings
        self.make_folder(settings)
        self.flaimapper(settings,first_task)
        self.test_errors(settings)

    def make_folder(self,settings):
        '''
        Make output folder for task if needed
        '''
        if not os.path.exists(os.path.join(settings["--output"],"identify")):
            os.makedirs(os.path.join(settings["--output"],"identify"))
        if not os.path.exists(os.path.join(settings["--output"],"identify","bam")):
            os.makedirs(os.path.join(settings["--output"],"identify","bam"))
        if not os.path.exists(os.path.join(settings["--output"],"identify","identify_info")):
            os.makedirs(os.path.join(settings["--output"],"identify","identify_info"))
        if not os.path.exists(os.path.join(settings["--output"],"identify","flaimapper")):
            os.makedirs(os.path.join(settings["--output"],"identify","flaimapper"))
        if not os.path.exists(os.path.join(settings["--output"],"identify","flaimapper","flaimapper_temp")):
            os.makedirs(os.path.join(settings["--output"],"identify","flaimapper","flaimapper_temp"))
        if not os.path.exists(os.path.join(settings["--output"],"identify","flaimapper","flaimapper_info")):
            os.makedirs(os.path.join(settings["--output"],"identify","flaimapper","flaimapper_info"))
        if not os.path.exists(os.path.join(settings["--output"],"identify","featurecounts")):
            os.makedirs(os.path.join(settings["--output"],"identify","featurecounts"))

    def flaimapper(self,settings,first_task):
        '''
        Make bam files and fragment them.
        Use Flaimapper on fragmented bam.
        Combine Flaimapper output and count reads using Bedtools intersect.
        '''
        overlap = settings["overlap_range"]
        size_range = settings["size_range"]
        #get chrom lengths
        genome = Fasta(settings["genome"],one_based_attributes=False)
        genome_lengths = {}
        for chrom in genome.keys():
            genome_lengths[chrom] = len(genome[chrom])
        #This allowes 2bp non-overlaping pp till length 1418
##        overlap = [0.888,0.916,0.939,0.957,0.97,0.979,0.985,0.9898,0.9931,\
##                   0.9953,0.99685,0.99789] #as list [x,y]
##        size_range = [24,33,47,67,97,140,197,292,432,636,950]\
##                     #as list [z], has to be shorter by 1 from overlap
        pool = mp.Pool(processes=settings["CPUs"])
        results = [pool.apply_async(self.flaimapper_by_library, \
                        args = (settings,library,overlap,size_range,\
                                genome_lengths,first_task)) \
                        for library in sorted(settings["libraries"])]
        pool.close()
        pool.join()
        for r in results:
            r.get()

    def flaimapper_by_library(self,settings,library,overlap,size_range,\
                              genome_lengths,first_task):
        '''
        Flaimapper by library.
        '''
        #make folders
        if not os.path.exists(os.path.join(settings["--output"],"identify","flaimapper",\
                                           "flaimapper_temp",library)):
            os.makedirs(os.path.join(settings["--output"],"identify","flaimapper",\
                                     "flaimapper_temp",library))
        if not os.path.exists(os.path.join(settings["--output"],"identify","flaimapper",\
                                           "flaimapper_info",library)):
            os.makedirs(os.path.join(settings["--output"],"identify","flaimapper",\
                                     "flaimapper_info",library))
            
        #set filenames
        if first_task == "identify":
            input_file = os.path.join(settings["--input"],library+\
                                        settings["identify"]["identify_input_file_suffix"])
        else:
            input_file = os.path.join(settings["--output"],"pseudoSE",\
                                          library+"_pseudoSE.sam")

        #set parameters
        strand_list = ["For", "Rev"]
        ##parameters to split reads by length
##        min_lower_range_limit = max(round(settings["min_length"],-1)-\
##                                settings["identify"]["identify_split_step"],0)
##        max_lower_range_limit = round(settings["max_length"],-1)-\
##                                settings["identify"]["identify_split_step"]
##        min_upper_range_limit = max(round(settings["min_length"],-1)+5,0)
##        max_upper_range_limit = round(settings["max_length"],-1)+5
##        length_range_upper = list(range(min_upper_range_limit,
##                                        max_upper_range_limit,\
##                                settings["identify"]["identify_split_step"]))
##        length_range_lower = list(range(min_lower_range_limit,
##                                        max_lower_range_limit,\
##                                settings["identify"]["identify_split_step"]))
        min_lower_range_limit = round(settings["min_length"],-1)
        if min_lower_range_limit > settings["min_length"]:
            min_lower_range_limit -= 10
        max_lower_range_limit = round(settings["max_length"],-1)
        min_upper_range_limit = min_lower_range_limit+settings["identify"]["identify_split_step"]+5
        max_upper_range_limit = max_lower_range_limit+5
        length_range_upper = list(range(min_upper_range_limit,
                                        max_upper_range_limit+settings["identify"]["identify_split_step"],\
                                settings["identify"]["identify_split_step"]))
        length_range_lower = list(range(min_lower_range_limit,
                                        max_lower_range_limit,\
                                settings["identify"]["identify_split_step"]))

        #split sam by parsing
        ##create files
        f_out_files_for = []
        f_out_files_rev = []
        for i in range(len(length_range_upper)):
            f_out_files_for.append(open(os.path.join(settings["--output"],"identify","bam",\
                                             library+"_"+"For"+"_"+str(i)+"_unsorted.sam"),"w"))
            f_out_full_files_for = open(os.path.join(settings["--output"],"identify","bam",\
                                             library+"_"+"For"+"_unsorted.sam"),"w")
            f_out_files_rev.append(open(os.path.join(settings["--output"],"identify","bam",\
                                             library+"_"+"Rev"+"_"+str(i)+"_unsorted.sam"),"w"))
            f_out_full_files_rev = open(os.path.join(settings["--output"],"identify","bam",\
                                             library+"_"+"Rev"+"_unsorted.sam"),"w")
            
        with open(input_file) as f_in:
            for line in f_in:
                #write header
                while line[0] == "@":
                    for i in range(len(length_range_upper)):
                        f_out_files_rev[i].write(line)
                        f_out_files_for[i].write(line)
                    f_out_full_files_for.write(line)
                    f_out_full_files_rev.write(line)
                    line = f_in.readline()

                length = int(line.split("\t")[5].strip("M"))
                index1 = bisect.bisect_left(length_range_lower,length)
                index2 = bisect.bisect_left(length_range_upper,length)
                
                if int(line.split("\t")[1]) & 16: #reverse
                    for i in range(index2,index1):
                        f_out_files_rev[i].write(line)
                    f_out_full_files_rev.write(line)
                else:
                    for i in range(index2,index1):
                        f_out_files_for[i].write(line)
                    f_out_full_files_for.write(line)
                        
        for i in range(len(length_range_upper)):
            f_out_files_rev[i].close()
            f_out_files_for[i].close()
        f_out_full_files_rev.close()
        f_out_full_files_for.close()

        #sam to bam
        for strand_name in strand_list:
            for i in range(len(length_range_upper)):    
                #sam to bam
                samtools_sam_to_bam = (
                    settings["samtools_call"], "view",
                    "-h", "-Sb", os.path.join(settings["--output"],"identify","bam",\
                                             library+"_"+strand_name+"_"+str(i)+"_unsorted.sam"),
                    "-o",os.path.join(settings["--output"],"identify","bam",\
                                             library+"_"+strand_name+"_"+str(i)+"_unsorted.bam")
                    )
                os.system(" ".join(samtools_sam_to_bam))

                #sort bam
                length_split_bam = os.path.join(settings["--output"],"identify","bam",\
                                             library+"_"+strand_name+"_"+str(i)+".bam")
                samtools_sort_command = (
                    settings["samtools_call"], "sort",
                    "-o", length_split_bam,
                    os.path.join(settings["--output"],"identify","bam",\
                                             library+"_"+strand_name+"_"+str(i)+"_unsorted.bam")
                    )

                os.system("\t".join(samtools_sort_command))
                                  
                #index bam
                samtools_index_command = (
                    settings["samtools_call"], "index",
                    length_split_bam
                    )

                os.system(" ".join(samtools_index_command))

                #remove unsorted bam
                os.remove(os.path.join(settings["--output"],"identify","bam",\
                                             library+"_"+strand_name+"_"+str(i)+"_unsorted.sam"))
                os.remove(os.path.join(settings["--output"],"identify","bam",\
                                             library+"_"+strand_name+"_"+str(i)+"_unsorted.bam"))
                
                #flaimapper
                flaimapper_output = os.path.join(settings["--output"],"identify","flaimapper",\
                                                 "flaimapper_temp",library,\
                                               library+"_"+strand_name+"_"+str(i)+"_flaimapper.tab")
                flaimapper_info = os.path.join(settings["--output"],"identify","flaimapper",\
                                               "flaimapper_info",library,\
                                               library+"_"+strand_name+"_"+str(i)+"_flaimapper.info")
                
                flaimapper_command = (
                    settings["identify"]["identify_flaimapper_call"],"-q", "-f 1",
                    "-o", flaimapper_output, "-q",
                    "-r", settings["genome"],
                    "-p", settings["identify"]["identify_flaimapper_parameters"],
                    length_split_bam, "2>" ,flaimapper_info
                    )
                #print(" ".join(flaimapper_command))
                os.system(" ".join(flaimapper_command))
                                  
                #remove length filtered bam and its indexes
                os.remove(length_split_bam)
                os.remove(length_split_bam+".bai")

            #format full strand files
            #sam to bam
            samtools_sam_to_bam = (
                settings["samtools_call"], "view",
                "-h", "-Sb", os.path.join(settings["--output"],"identify","bam",\
                                         library+"_"+strand_name+"_unsorted.sam"),
                "-o",os.path.join(settings["--output"],"identify","bam",\
                                         library+"_"+strand_name+"_unsorted.bam")
                )
            os.system(" ".join(samtools_sam_to_bam))

            #sort bam
            strand_split_bam = os.path.join(settings["--output"],"identify","bam",\
                                         library+"_"+strand_name+".bam")
            samtools_sort_command = (
                settings["samtools_call"], "sort",
                "-o", strand_split_bam,
                os.path.join(settings["--output"],"identify","bam",\
                                         library+"_"+strand_name+"_unsorted.bam")
                )

            os.system("\t".join(samtools_sort_command))
                              
            #index bam
            samtools_index_command = (
                settings["samtools_call"], "index",
                strand_split_bam
                )

            os.system(" ".join(samtools_index_command))

            #remove unsorted bam
            os.remove(os.path.join(settings["--output"],"identify","bam",\
                                         library+"_"+strand_name+"_unsorted.sam"))
            os.remove(os.path.join(settings["--output"],"identify","bam",\
                                         library+"_"+strand_name+"_unsorted.bam"))


        
##        #process by strand
##        strand_split_parameter = ["-F","-f"] #needed for spliting by strand
##        for j,strand_name in enumerate(strand_list):
##            print("\t"+" ".join([library,strand_name,"Flaimapper"]))
##            
##            #split by strand
##            strand_split_unsorted_bam = os.path.join(settings["--output"],"identify",\
##                                             "bam",library+"_"+strand_name+"_unsorted.bam")
##            samtools_split_by_strand_command = (
##                        settings["samtools_call"], "view",
##			"-@", str(settings["samtools_threads"]),
##                        strand_split_parameter[j], "0x10 -Sb",
##                        input_file, ">",
##                        strand_split_unsorted_bam
##                        )
##            os.system("\t".join(samtools_split_by_strand_command))
##            
##            #sort stranded bam
##            strand_split_bam = os.path.join(settings["--output"],"identify","bam",\
##                                             library+"_"+strand_name+".bam")
##            samtools_sort_command = (
##                settings["samtools_call"], "sort",
###				"-@", str(settings["samtools_threads"]),
##                "-o", strand_split_bam, strand_split_unsorted_bam
##                )
##            os.system("\t".join(samtools_sort_command))
##
##            #index bam
##            samtools_index_command = (
##                settings["samtools_call"], "index",
##				"-@", str(settings["samtools_threads"]),
##                strand_split_bam
##                )
##
##            os.system(" ".join(samtools_index_command))
##            
##            #sorted bam to sam
##            samtools_bam_to_sam_command = (
##                settings["samtools_call"], "view", "-h",
##				"-@", str(settings["samtools_threads"]),
##                "-o", strand_split_bam[:-3]+"sam", strand_split_bam
##                )
##            os.system("\t".join(samtools_bam_to_sam_command))              
##            #remove unsorted stranded bam
##            os.remove(strand_split_unsorted_bam)
##            
##            #split bam by length
##            for i in range(len(length_range_upper)):
##                length_split_unsorted_bam = os.path.join(settings["--output"],"identify",\
##                                                 "bam",library+\
##                                               "_"+strand_name+"_"+str(i)+"_unsorted.bam")
##                samtools_split_by_length_command = (
##                    settings["samtools_call"], "view",
##					"-@", str(settings["samtools_threads"]),
##                    "-h", strand_split_bam,
##                    "| awk 'length($10) >",
##                    str(length_range_lower[i]),
##                    "&& length($10) <",
##                    str(length_range_upper[i]),
##                    "|| $1 ~ /^@/' |",
##                    settings["samtools_call"], "view",
##                    "-bS - >", length_split_unsorted_bam
##                    )
##                #print(" ".join(samtools_split_by_length_command))
##                os.system(" ".join(samtools_split_by_length_command))
##
##
##                #sort bam
##                length_split_bam = os.path.join(settings["--output"],"identify","bam",\
##                                                 library+"_"+strand_name+"_"+str(i)+".bam")
##                samtools_sort_command = (
##                    settings["samtools_call"], "sort",
##					"-@", str(settings["samtools_threads"]),
##                    "-o", length_split_bam, length_split_unsorted_bam
##                    )
##
##                os.system("\t".join(samtools_sort_command))
##                                  
##                #index bam
##                samtools_index_command = (
##                    settings["samtools_call"], "index",
##					"-@", str(settings["samtools_threads"]),
##                    length_split_bam
##                    )
##
##                os.system(" ".join(samtools_index_command))
##
##                #remove unsorted bam
##                os.remove(length_split_unsorted_bam)
##                
##                #flaimapper
##                flaimapper_output = os.path.join(settings["--output"],"identify","flaimapper",\
##                                                 "flaimapper_temp",library,\
##                                               library+"_"+strand_name+"_"+str(i)+"_flaimapper.tab")
##                flaimapper_info = os.path.join(settings["--output"],"identify","flaimapper",\
##                                               "flaimapper_info",library,\
##                                               library+"_"+strand_name+"_"+str(i)+"_flaimapper.info")
##                
##                flaimapper_command = (
##                    settings["identify"]["identify_flaimapper_call"],"-q", "-f 1",
##                    "-o", flaimapper_output, "-q",
##                    "-r", settings["genome"],
##                    "-p", settings["identify"]["identify_flaimapper_parameters"],
##                    length_split_bam, "2>" ,flaimapper_info
##                    )
##                #print(" ".join(flaimapper_command))
##                os.system(" ".join(flaimapper_command))
##                                  
##                #remove length filtered bam and its indexes
##                os.remove(length_split_bam)
##                os.remove(length_split_bam+".bai")

            #COMBINE FLAIMAPPER OUTPUT
            ##READ IN DATA
            list_pp,stat_pp = self.read_in_flaimapper_data(settings,\
                                    len(length_range_upper),strand_name,library)
                      
            ##WRITE DATA INTO FILE
            print("\t"+" ".join([library,"Write data"]))
            file_name = os.path.join(settings["--output"],"identify",\
                                     library+"_"+strand_name+"_pp")
            infofile_name = os.path.join(settings["--output"],"identify",\
                                    "identify_info",\
                                     library+"_"+strand_name+"_identifyinfo.log")
            self.write_BED(file_name+".BED",list_pp,genome_lengths)
            self.write_SAF(file_name+".SAF",list_pp,genome_lengths)
            self.write_statistics(infofile_name,list_pp,stat_pp)
            
            #COUNT READS PER PP
            print("\t"+" ".join([library,"Count reads"]))
##            self.fragment_BED(file_name+".BED",overlap,size_range)
##            self.count_reads_fragmented_BED_bedtoools(\
##                settings,file_name+".BED",strand_split_bam,library,overlap,strand_name)
##            self.remove_fragmented_BED(file_name+".BED",overlap)
            ##count reads by featurecounts
            #self.fragment_SAF(file_name+".SAF",overlap,size_range)
            #self.count_reads_fragmented_SAF_featurecounts(\
            #    settings,file_name+".SAF",strand_split_bam,library,overlap,strand_name)
            self.count_reads_featurecounts(\
                settings,file_name+".SAF",strand_split_bam,library,overlap,strand_name)
            #self.remove_fragmented_SAF(file_name+".SAF",overlap)
        
    def read_in_flaimapper_data(self,settings,input_file_range,strand_name,library):
        '''
        Reads in and combines Flaimapper data
        '''
        #set parameters
        set_pp = set()
        stat_pp = {}
        stat_pp["length"] = {}
        stat_pp["count"] = {}
        if strand_name == "For":
            strand = "+"
        elif strand_name == "Rev":
            strand = "-"
        #f_out_counts = open(output_flaimapper_counts,"w") #for observing flaimapper pp-s
        print("\t"+" ".join([library,strand_name,"Read in data"]))
        for i in range(input_file_range):
            input_file = os.path.join(settings["--output"],"identify","flaimapper",\
                                      "flaimapper_temp",library,\
                                        library+"_"+strand_name+"_"+str(i)+"_flaimapper.tab")
                                      
            #The first base in a chromosome is numbered 0 in flaimapper tabular output,
            #end inclusive
            f_in = open(input_file)
            f_in.readline() #to remove description line
            for line in f_in:
                line = line.strip().split("\t")
                chrom = line[2]
                start = int(line[3])
                end = int(line[4])
                #flaimapper gives reads for both ends, take smallest
                read_number = min(int(line[9]),int(line[10])) 
                length = int(line[1])
                #excluding flaimapper pp-s by various parameters
                if length < settings["min_length"]:
                    continue
                if length > settings["max_length"]:
                    continue
                if read_number < settings["min_pp_reads"]/2:
                    continue
                #collect statistics only when pp not represented allready
                if (chrom,start,end,strand) not in set_pp:
                    if read_number not in stat_pp["count"]:
                        stat_pp["count"][read_number] = 1
                    else:
                        stat_pp["count"][read_number] += 1
                    if length not in stat_pp["length"]:
                        stat_pp["length"][length] = 1
                    else:
                        stat_pp["length"][length] += 1
                set_pp.add((chrom,start,end,strand))
            f_in.close()
        list_pp = sorted(list(set_pp))
        return list_pp,stat_pp

    def write_BED(self,output_BED,list_pp,genome_lengths):
        '''
        Writes combined and selected flaimapper output into BED format.
        '''
        #write to BED
        #The first base in a chromosome is numbered 0 (0-based)
        #end exclusive
        f_out_BED = open(output_BED,"w")
        f_out_BED.write("#Chr\tStart\tEnd\tGeneID\tScore\tStrand\n")
        if len(list_pp) > 0:
            #last_pp_start = list_pp[-1][1]
            for pp in list_pp:
                f_out_BED.write(pp[0]+"\t")
                f_out_BED.write(str(pp[1])+"\t")
                f_out_BED.write(str(pp[2]+1)+"\t") #end exclusive
                f_out_BED.write("FM_"+pp[3]+pp[0]+"_"+(len(str(genome_lengths[pp[0]]))-\
                                len(str(pp[1])))*"0"+str(pp[1])+"_"+str(pp[2]-pp[1]+1)+"\t")
                f_out_BED.write(str(0)+"\t")
                f_out_BED.write(pp[3]+"\n")
        f_out_BED.close()

    def write_SAF(self,output_SAF,list_pp,genome_lengths):
        #write in to SAF
        #The first base in a chromosome is numbered 1 (1-based)
        #end inclusive
        f_out_pp = open(output_SAF,"w")
        f_out_pp.write("#GeneID\tChr\tStart\tEnd\tStrand\n")
        if len(list_pp) > 0:
            last_pp_start = list_pp[-1][1]
            for pp in list_pp:
                f_out_pp.write("FM_"+pp[3]+pp[0]+"_"+(len(str(genome_lengths[pp[0]]))-\
                               len(str(pp[1])))*"0"+str(pp[1])+"_"+str(pp[2]-pp[1]+1)+"\t")
                f_out_pp.write(pp[0]+"\t")
                f_out_pp.write(str(pp[1]+1)+"\t")
                f_out_pp.write(str(pp[2]+1)+"\t")
                f_out_pp.write(pp[3]+"\n")
        f_out_pp.close()
            
    def write_statistics(self,output_stat,list_pp,stat_pp):
        '''
        Writes statistics about combined and selected flaimapper output.
        '''
        f_out_stat = open(output_stat,"w")
        f_out_stat.write("Total "+str(sum([stat_pp["length"][x] for x in stat_pp["length"]])) + " processing products detected\n")
        f_out_stat.write("\n")
        f_out_stat.write("Length distribution:\n")
        for length in sorted(stat_pp["length"]):
            f_out_stat.write("\t"+str(length)+":\t"+str(stat_pp["length"][length])+"\n")

        f_out_stat.write("\n")
        f_out_stat.write("Count distribution:\n")
        for count in sorted(stat_pp["count"]):
            f_out_stat.write("\t"+str(count)+":\t"+str(stat_pp["count"][count])+"\n")
        f_out_stat.close()                            
                                                   
    def fragment_BED(self,input_BED,overlap,size_range):
        '''
        Fragments input_BED file by size ranges given
        '''
        f_BED_in = open(input_BED)
        #create temporary output files
        f_out = []
        f_out_names = []
        for i in range(len(overlap)):
            f_out.append(open(input_BED[:-4]+"_"+str(i)+".BED","w"))
            f_out_names.append(input_BED[:-4]+"_"+str(i)+".BED")    
        #go through input file and write pp-s to temporary files
        for line in f_BED_in:
            if line[0] == "#": #skip header
                line = f_BED_in.readline()
                if len(line) == 0: #stop if file is empty
                    break
            split_line = line.strip().split("\t")
            #gets index for output file
            #in short it will compares length of the pp with the size range and gives index for that
            index = bisect.bisect(size_range, int(split_line[2])-int(split_line[1]))             
            f_out[index].write(line)
        #closing output files
        for i in range(len(overlap)):
            f_out[i].close()
        f_BED_in.close()
        
    def count_reads_fragmented_BED_bedtoools(self,settings,input_BED,input_bam,library,overlap,strand_name):
        '''
        Count reads for fragmented input_BED file by bestools intersect.
        COmbine result files afterwords
        '''
                                         
        #count reads in different files           
        for i in range(len(overlap)):
            count_reads_command =(
                    settings["bedtools_call"], "intersect",
                    "-a", input_BED[:-4]+"_"+str(i)+".BED",\
                    "-b", input_bam,
                    "-f", str(overlap[i]),\
                    "-r", "-sorted","-s", "-c",\
                    ">", input_BED[:-4]+"_"+str(i)+"_counted.BED"
                    )
            os.system(" ".join(count_reads_command))
         
        #combine fragmented_pp_counted_files 
        with open(input_BED[:-4]+"_counted_unsorted.BED",'wb') as wfd:
            for i in range(len(overlap)):
                with open(input_BED[:-4]+"_"+str(i)+"_counted.BED",'rb') as fd:
                    shutil.copyfileobj(fd, wfd, 1024*1024*10)
                        
        #sort combined counted pp-s 
        sort_command = (
                settings["bedtools_call"], "sort",\
                "-i",input_BED[:-4]+"_counted_unsorted.BED",\
                ">",input_BED[:-4]+"_counted.BED"
                )
        os.system(" ".join(sort_command))
            
        #delete temporary count files         
        for i in range(len(overlap)):
            os.remove(input_BED[:-4]+"_"+str(i)+"_counted.BED") #remove fragmented_pp_counted_files
        os.remove(input_BED[:-4]+"_counted_unsorted.BED") #remove unsorted combined file
        
    def remove_fragmented_BED(self,input_BED,overlap):
        '''
        Remove fragmented input_BED files
        '''
        #delete fragmented input_BED_file parts       
        for i in range(len(overlap)):
            os.remove(input_BED[:-4]+"_"+str(i)+".BED") #remove fragmented_pp_counted_files

    def fragment_SAF(self,input_SAF,overlap,size_range):
        '''
        Fragments input_SAF file by size ranges given
        '''
        ##fragment input_BED_file by size range given
        f_SAF_in = open(input_SAF)
        #create temporary output files
        f_out = []
        f_out_names = []
        for i in range(len(overlap)):
            f_out.append(open(input_SAF[:-4]+"_"+str(i)+".SAF","w"))
            f_out_names.append(input_SAF[:-4]+"_"+str(i)+".SAF")    
        #go through input file and write pp-s to temporary files
        for line in f_SAF_in:
            if line[0] == "#": #skip header
                line = f_SAF_in.readline()
                if len(line) == 0: #stop if file is empty
                    break
            split_line = line.strip().split("\t")
            #gets index for output file
            #in short it will compares length of the pp with the size range and gives index for that
            index = bisect.bisect(size_range, int(split_line[3])-int(split_line[2])+1)             
            f_out[index].write(line)
        #closing output files
        for i in range(len(overlap)):
            f_out[i].close()
        f_SAF_in.close()

    def count_reads_fragmented_SAF_featurecounts(self,settings,input_SAF,input_bam,library,overlap,strand_name):
        '''
        Count reads for fragmented input_BED file by featureCounts
        '''
        #count reads in different files
        for i in range(len(overlap)):
            if not os.path.exists(os.path.join(settings["--output"],"identify","flaimapper",str(i))):
                os.makedirs(os.path.join(settings["--output"],"identify","flaimapper",str(i)))
            if os.stat(input_SAF[:-4]+"_"+str(i)+".SAF").st_size != 0:
                feturecounts_info = os.path.join(settings["--output"],"identify","featurecounts",\
                                   library+"_"+strand_name+"_"+str(i)+"_featurecounts.info")
                file_name = os.path.join(settings["--output"],"identify","flaimapper",str(i),\
                                library+"_"+strand_name+"_pp") 
                featureCounts_command =(
                                settings["featureCounts_call"],
                                "-T", str(settings["CPUs"]),
                                #"-G", settings["genome"],
                                "-M", "-O", "-s 1", "-F SAF", "-R SAM",
                                "--fracOverlap", str(overlap[i]),
                                "--fracOverlapFeature", str(overlap[i]),                                             
                                "-a", input_SAF[:-4]+"_"+str(i)+".SAF",
                                "-o", file_name+"_"+str(i)+"_counted.SAF",
                                input_bam, "2>", feturecounts_info
                                )       
                os.system(" ".join(featureCounts_command))
         
        #combine fragmented_pp_counted_files 
        with open(input_SAF[:-4]+"_counted_unsorted.SAF",'wb') as wfd:
            for i in range(len(overlap)):
                if os.stat(input_SAF[:-4]+"_"+str(i)+".SAF").st_size != 0:
                    with open(input_SAF[:-4]+"_"+str(i)+"_counted.SAF",'rb') as fd:
                        shutil.copyfileobj(fd, wfd, 1024*1024*10)
                        
        #sort combined counted pp-s 
        sort_command = (
                settings["bedtools_call"], "sort",\
                "-i",input_SAF[:-4]+"_counted_unsorted.SAF",\
                ">",input_SAF[:-4]+"_counted.SAF"
                )
        os.system(" ".join(sort_command))
            
        #delete temporary count files         
        #for i in range(len(overlap)):
        #    os.remove(input_SAF[:-4]+"_"+str(i)+"_counted.SAF") #remove fragmented_pp_counted_files
        #os.remove(input_SAF[:-4]+"_counted_unsorted.SAF") #remove unsorted combined file         

    def count_reads_featurecounts(self,settings,input_SAF,input_bam,library,overlap,strand_name):
        '''
        Count reads input_SAF file by featureCounts
        '''
        #count reads in different files
        featurecounts_info = os.path.join(settings["--output"],"identify","featurecounts",\
                           library+"_"+strand_name+"_featurecounts.info")
        output_SAF =input_SAF[:-4]+"_counted.SAF"

        featureCounts_command =(
                        settings["featureCounts_call"],
                        #"-T", str(settings["CPUs"]),
                        #"-G", settings["genome"],
                        "-M", "-O", "-s 1", "-F SAF",
                        "--nonOverlap", str(settings["non_overlap"]),
                        "--nonOverlapFeature", str(settings["non_overlap"]),
##                        "--fracOverlap", str(overlap[i]),
##                        "--fracOverlapFeature", str(overlap[i]),                                             
                        "-a", input_SAF,
                        "-o", output_SAF,
                        input_bam, "2>", featurecounts_info
                        )
        os.system(" ".join(featureCounts_command))
        if os.path.isfile(output_SAF):
            print("AA",output_SAF)
        else:
            print("BB",output_SAF)        
##        #combine fragmented_pp_counted_files 
##        with open(input_SAF[:-4]+"_counted_unsorted.SAF",'wb') as wfd:
##            for i in range(len(overlap)):
##                if os.stat(input_SAF[:-4]+"_"+str(i)+".SAF").st_size != 0:
##                    with open(input_SAF[:-4]+"_"+str(i)+"_counted.SAF",'rb') as fd:
##                        shutil.copyfileobj(fd, wfd, 1024*1024*10)
        
        #convert SAF to BED
        #this is done in historical reasons just not to change cluster.py and
        #decrease number of file formats used
        self.SAF_to_BED(output_SAF,input_SAF[:-4]+"_counted.BED")
             
##        #sort combined counted pp-s
##        sort_command = (
##                settings["bedtools_call"], "sort",\
##                "-i",input_SAF[:-4]+"_counted_unsorted.BED",\
##                ">",input_SAF[:-4]+"_counted.SAF"
##                )
####        sort_command = (
####            "sort", "-k2,2", "-k3n,3", "-k4n,4",
####            "-t", "\t",
####            input_SAF[:-4]+"_counted_unsorted.SAF",
####            ">",input_SAF[:-4]+"_counted.SAF"
####            )
##        os.system(" ".join(sort_command))
            
        #delete temporary count files         
        #for i in range(len(overlap)):
        #    os.remove(input_SAF[:-4]+"_"+str(i)+"_counted.SAF") #remove fragmented_pp_counted_files
        #remove SAF files
        #os.remove(output_SAF) #remove counted SAF file
        #os.remove(input_SAF) #remove SAF file
        #os.remove(input_SAF[:-4]+"_counted_unsorted.SAF.summary") #remove unsorted combined file
        #os.remove(input_SAF[:-4]+"_counted_unsorted.BED") #remove unsorted combined file 

    def SAF_to_BED(self,SAF_file,BED_file):
        '''
        Converts SAF to BED.
        '''
        if os.path.isfile(SAF_file):
            print("A",SAF_file)
        else:
            print("B",SAF_file)
        f_out = open(BED_file,"w")
        print(sorted(os.listdir(os.path.dirname(SAF_file))))
        with open(SAF_file) as f_in:
            #skip featureCounts command and header lines
            f_in.readline()
            f_in.readline()
            for line in f_in:
                line = line.strip().split("\t")
                f_out.write("\t".join([line[1],str(int(line[2])-1),line[3],line[0],"0",line[4],line[6]])+"\n")
        f_out.close()
        
    def remove_fragmented_SAF(self,input_SAF,overlap):
        '''
        Remove fragmented input_BED files
        '''
        #delete fragmented input_BED_file parts       
        #for i in range(len(overlap)):
        #    os.remove(input_SAF[:-4]+"_"+str(i)+".SAF") #remove fragmented_pp_counted_files


    def test_errors(self,settings):
        '''
        Testing existence of output files
        '''
        output_folder = os.path.join(settings["--output"],"identify")
        files = os.listdir(output_folder)
        if len(files) == 0:
            sys.exit("'Task' "+pseudoSE+" incomplete!")
        for library in settings["libraries"]:
            for suffix in ["_For_pp.BED","_For_pp_counted.BED","_Rev_pp.BED","_Rev_pp_counted.BED"]:
                if not (library + suffix) in files:
                    sys.exit('Error in task "identify", problem with library: ' + library)

