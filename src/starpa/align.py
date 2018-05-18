#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
To align reads to genome with bowtie2

To do:

'''

import os
import sys
import multiprocessing as mp

class align():
    def __init__(self,settings,mode,first_task):
        self.settings = settings
        self.make_folder(settings,mode)
        self.make_index(settings)
        if settings["paired"]:
            pool = mp.Pool(processes=settings["CPUs"])
            results = [pool.apply_async(self.bowtie2_PE, \
                              args = (settings,library,mode,first_task)) \
                              for library in sorted(settings["libraries"])]
            pool.close()
            pool.join()
            for r in results:
                r.get()
        else:
            pool = mp.Pool(processes=settings["CPUs"])
            results = [pool.apply_async(self.bowtie2_SE, \
                              args = (settings,library,mode,first_task)) \
                              for library in sorted(settings["libraries"])]
            pool.close()
            pool.join()
            for r in results:
                r.get()
        if mode == "II":
            self.combine_log_files(settings)
        self.test_errors(settings,mode)

    def make_folder(self,settings,mode):
        '''
        Make output folder for task if needed
        '''
        if not os.path.exists(os.path.join(settings["--output"],"align")):
            os.makedirs(os.path.join(settings["--output"],"align"))
        if mode != "all":
            if not os.path.exists(os.path.join(settings["--output"],"align","round_"+mode)):
                os.makedirs(os.path.join(settings["--output"],"align","round_"+mode))
            if not os.path.exists(os.path.join(settings["--output"],"align","round_"+\
                                               mode,"align_info")):
                os.makedirs(os.path.join(settings["--output"],"align","round_"+\
                                         mode,"align_info"))
        if mode != "I":
            if not os.path.exists(os.path.join(settings["--output"],"align","align_info")):
                os.makedirs(os.path.join(settings["--output"],"align","align_info"))

    def make_index(self,settings):
        '''
        Run bowtie2-build to make index of reference sequence
        '''
        #test existence of bowtie2 index files
        index_files = ["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]       
        files = os.listdir(os.path.join("/".join(settings["genome"].split("/")[:-1])))
        for index_file in index_files:
            if not any(file.endswith(index_file) for file in files):
                bowtie2_index_command = (
                            settings["align"]["align_index_call"],
                            "-f", settings["genome"], "-q",
                            ".".join(settings["genome"].split(".")[:-1])
                            )
                os.system(" ".join(list(bowtie2_index_command)))
                break
            
    def bowtie2_PE(self,settings,library,mode,first_task):
        '''
        Run bowtie2 with paired end settings
        '''
        print("\t"+library)
        #set filenames
        if mode != "all":
            if mode == "II":
                input_file1 = os.path.join(settings["--output"],"align","round_I",\
                                          library+"_unmapped_1.fq")
                input_file2 = os.path.join(settings["--output"],"align","round_I",\
                                          library+"_unmapped_2.fq")
            else:
                if first_task == "align":
                    input_file1 = os.path.join(settings["--input"],library+\
                                               settings["align"]["align_input_file_suffix_for"])
                    input_file2 = os.path.join(settings["--input"],library+\
                                               settings["align"]["align_input_file_suffix_rev"])
                else:
                    input_file1 = os.path.join(settings["--output"],"trim",library+\
                                               settings["align"]["align_input_file_suffix_for"])
                    input_file2 = os.path.join(settings["--output"],"trim",library+\
                                               settings["align"]["align_input_file_suffix_rev"])
                
            output_file = os.path.join(settings["--output"],"align","round_"+mode,library+\
                                       "_"+mode+".sam")
            align_info_file = os.path.join(settings["--output"],"align","round_"+mode,\
                                      "align_info",library+"_"+mode+"_aligninfo.log")       

        else:
            if first_task == "align":
                input_file1 = os.path.join(settings["--input"],library+\
                                           settings["align"]["align_input_file_suffix_for"])
                input_file2 = os.path.join(settings["--input"],library+\
                                           settings["align"]["align_input_file_suffix_rev"])                    
            else:
                input_file1 = os.path.join(settings["--output"],"trim",library+\
                                           settings["align"]["align_input_file_suffix_for"])
                input_file2 = os.path.join(settings["--output"],"trim",library+\
                                           settings["align"]["align_input_file_suffix_rev"])
            output_file = os.path.join(settings["--output"],"align",library+".sam")
            align_info_file = os.path.join(settings["--output"],"align",\
                                      "align_info",library+"_aligninfo.log")
       
        #set phred score type
        if str(settings["trim"]["trim_quality_base"]) == "64":
            phred_type = "--phred64"
        else:
            phred_type = "--phred33"
        #bowtie2 out put is phred 33
        if mode == "II":
            phred_type = "--phred33"
            
        #set -L parameter
        if mode == "II":
            seed_length = "14"
        else:
            seed_length = "22"
            
        bowtie2_command = (
                    settings["align"]["align_call"],
                    "-x", ".".join(settings["genome"].split(".")[:-1]), #removes extention
                                                                        #from file name
                    "-p", str(settings["align"]["align_threads"]),
                    "--rdg", "100,3",  #to avoid indels
                    "--rfg", "100,3",  #to avoid indels
                    phred_type,
                    "-N", "1", #increased sensitivity
                    "-a", #unlimited number of matches
                    "-L", seed_length,
                    "--no-mixed", "--no-discordant",
                    "-X", str(settings["max_length"]),
                    "-I", str(settings["min_length"]),
                    "--score-min","L,0,-2", #alignments with more mismatches are allowed
                    "--rg-id", "N1_"+phred_type.strip("--")+"_gap100",
                    "-1", input_file1, "-2", input_file2,
                    ">", output_file, "2>", align_info_file
                    )
        os.system(" ".join(list(bowtie2_command)))
    
    def bowtie2_SE(self,settings,library,mode,first_task):
        '''
        Run bowtie2 with single end settings
        '''
        print("\t"+library)
        #set filenames       
        if mode != "all":
            if mode == "II":
                input_file = os.path.join(settings["--output"],"align","round_I",\
                                          library+"_unmapped.fq")
            else:
                if first_task == "align":
                    input_file = os.path.join(settings["--input"],library+\
                                              settings["align"]["align_input_file_suffix_SE"])
                else:
                    input_file = os.path.join(settings["--output"],"trim",library+\
                                              settings["align"]["align_input_file_suffix_SE"])
    
            output_file = os.path.join(settings["--output"],"align","round_"+mode,library+\
                                       "_"+mode+".sam")
            align_info_file = os.path.join(settings["--output"],"align","round_"+mode,\
                                      "align_info",library+"_"+mode+"_aligninfo.log")       

        else:
            if first_task == "align":
                input_file = os.path.join(settings["--input"],library+\
                                          settings["align"]["align_input_file_suffix_SE"])
            else:
                input_file = os.path.join(settings["--output"],"trim",library+\
                                          settings["align"]["align_input_file_suffix_SE"])
            output_file = os.path.join(settings["--output"],"align",library+".sam")
            align_info_file = os.path.join(settings["--output"],"align",\
                                      "align_info",library+"_aligninfo.log")
                
        #set phred score type
        if str(settings["trim"]["trim_quality_base"]) == 64:
            phred_type = "--phred64"
        else:
            phred_type = "--phred33"
            
        #set -L parameter
        if mode == "second":
            seed_length = "14"
        else:
            seed_length = "22"
        bowtie2_command = (
                    settings["align"]["align_call"],
                    "-x", ".".join(settings["genome"].split(".")[:-1]), #removes extention
                                                                        #from file name
                    "-p", str(settings["align"]["align_threads"]),
                    "--rdg", "100,3",  #to avoid indels
                    "--rfg", "100,3",  #to avoid indels
                    phred_type,
                    "-N", "1", #increased sensitivity
                    "-a", #unlimited number of matches
                    "-L", seed_length,
                    "--score-min","L,0,-2", #alignments with more mismatches are allowed
                    "--rg-id", "N1_phred64_gap100",
                    "-U", input_file,
                    ">", output_file, "2>", align_info_file
                    )
        #print("AA",bowtie2_command)
        os.system(" ".join(list(bowtie2_command)))

    def combine_log_files(self,settings):
        '''
        Combines log files from round I and round III
        '''
        for library in sorted(settings["libraries"]):
            #create common info file
            report_file1 = os.path.join(settings["--output"],"align","round_I",\
                                           "align_info",library+"_I_aligninfo.log")
            report_file2 = os.path.join(settings["--output"],"align","round_II",\
                                           "align_info",library+"_II_aligninfo.log")
            #create parameters
            total_reads = 0
            group_reads = 0
            not_aligned_reads = 0
            mono_aligned_reads = 0
            multi_aligned_reads = 0
            
            ##read in from file I
            with open(report_file1) as f_in1:
                line = f_in1.readline()
                total_reads = int(line.split(" ")[0])
                if line.find("; of these:") != -1:
                    group_reads = int(f_in1.readline().strip().split(" ")[0])
                    f_in1.readline()
                    #not_aligned_reads = int(f_in1.readline().strip().split(" ")[0])
                    mono_aligned_reads = int(f_in1.readline().strip().split(" ")[0])
                    multi_aligned_reads = int(f_in1.readline().strip().split(" ")[0])
                
            ##read in from file II
            with open(report_file2) as f_in2:
                line = f_in2.readline()
                if line.find("; of these:") != -1:
                    f_in2.readline()
                    #total_reads = int(f_in2.readline().split(" ")[0])
                    #paired_reads = int(f_in2.readline().strip().split(" ")[0])
                    not_aligned_reads = int(f_in2.readline().strip().split(" ")[0])
                    mono_aligned_reads += int(f_in2.readline().strip().split(" ")[0])
                    multi_aligned_reads += int(f_in2.readline().strip().split(" ")[0])
                
            ##write new reportfile
            report_file = os.path.join(settings["--output"],"align",\
                                           "align_info",library+"_aligninfo.log")                
            with open(report_file, "w") as f_out:
                f_out.write(str(total_reads) + " reads; of these:\n")
                if settings["paired"]: 
                    f_out.write("  "+str(group_reads)+" ("+\
                                str(round(group_reads/total_reads*100,2))+\
                                "%) were paired; of these:\n")
                    f_out.write("    "+str(not_aligned_reads)+" ("+\
                                str(round(not_aligned_reads/total_reads*100,2))+\
                                "%) aligned concordantly 0 times\n")
                    f_out.write("    "+str(mono_aligned_reads)+" ("+\
                                str(round(mono_aligned_reads/total_reads*100,2))+\
                                "%) aligned concordantly exactly 1 time\n")
                    f_out.write("    "+str(multi_aligned_reads)+" ("+\
                                str(round(multi_aligned_reads/total_reads*100,2))+\
                                "%) aligned concordantly > 1 time\n")
                else:
                    f_out.write("  "+str(group_reads)+" ("+\
                                str(round(group_reads/total_reads*100,2))+\
                                "%) were unpaired; of these:\n")                    
                    f_out.write("    "+str(not_aligned_reads)+" ("+\
                                str(round(not_aligned_reads/total_reads*100,2))+\
                                "%) aligned 0 times\n")
                    f_out.write("    "+str(mono_aligned_reads)+" ("+\
                                str(round(mono_aligned_reads/total_reads*100,2))+\
                                "%) aligned exactly 1 time\n")
                    f_out.write("    "+str(multi_aligned_reads)+" ("+\
                                str(round(multi_aligned_reads/total_reads*100,2))+\
                                "%) aligned > 1 time\n")
                f_out.write(str(round(group_reads/total_reads*100,2))+\
                            "% overall alignment rate\n")                    

    def test_errors(self,settings,mode):
        '''
        Tests are any info files containing an error raport
        '''
        #info folder names
        if mode != "all": 
            info_folder = os.path.join(settings["--output"],"align","round_"+mode,\
                                      "align_info")       
        else:
            info_folder = os.path.join(settings["--output"],"align",\
                                      "align_info")
        if len(os.listdir(info_folder)) == 0:
            sys.exit('Task "align" incomplete! Infofiles missing in folder ' + info_folder)
        #test content of info files
        for file in os.listdir(info_folder):
            with open(os.path.join(info_folder,file)) as f_in:
                for line in f_in:
                    if line.startswith("Warning:") or line.startswith("Error"):
                        sys.exit('Error in task "align", check infofile: '+file)
            
        
        
 
