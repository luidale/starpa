#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Sorts aligned reads.

To do:

'''

import os
import sys
import multiprocessing as mp

class sam_sort():
    def __init__(self,settings,mode,first_task):
        self.settings = settings
        #self.mode = mode
        self.make_folder(settings,mode)
        if settings["paired"]:
            self.sam_sort_PE(settings,mode,first_task)
            if mode == "I":
                self.sam_unmapped_to_fq_PE(settings)  
            elif mode == "II":
                self.join_sam(settings)
                
        else:
            self.sam_sort_SE(settings,mode,first_task)
            if mode == "I":
                self.sam_unmapped_to_fq_SE(settings)
            elif mode == "II":
                self.join_sam(settings)
        self.test_errors(settings,mode)

        
    def make_folder(self,settings,mode):
        '''
        Make output folder for task if needed
        '''
        if not os.path.exists(os.path.join(settings["--output"],"sam_sort")):
            os.makedirs(os.path.join(settings["--output"],"sam_sort"))
        if mode != "all":
            if not os.path.exists(os.path.join(settings["--output"],"sam_sort","round_"+mode)):
                os.makedirs(os.path.join(settings["--output"],"sam_sort","round_"+mode))
            if not os.path.exists(os.path.join(settings["--output"],"sam_sort","round_"+\
                                               mode,"sam_sort_info")):
                os.makedirs(os.path.join(settings["--output"],"sam_sort","round_"+\
                                         mode,"sam_sort_info"))
        if mode != "I":
            if not os.path.exists(os.path.join(settings["--output"],"sam_sort","sam_sort_info")):
                os.makedirs(os.path.join(settings["--output"],"sam_sort","sam_sort_info"))                   

    def sam_sort_PE(self,settings,mode,first_task):
        '''
        For paired-end data.
        Sorts out unmapped and mappings with low mapping quality.
        Leaving only best mappings of each read-pair.
        In case primary mapping was eliminated the sam flag and sequence quality
        fields are corrected.
        '''
        pool = mp.Pool(processes=settings["CPUs"])
        results = [pool.apply_async(self.sam_sort_PE_library, \
                          args = (settings,library,mode,first_task)) \
                          for library in sorted(settings["libraries"])]
        pool.close()
        pool.join()
        for r in results:
            r.get()

    def sam_sort_SE(self,settings,mode,first_task):
        '''
        For single-end data.
        Sorts out unmapped and mappings with low mapping quality.
        Leaving only best mappings of each read.
        In case primary mapping was eliminated the sam flag and sequence quality
        fields are corrected.
        '''
        pool = mp.Pool(processes=settings["CPUs"])
        results = [pool.apply_async(self.sam_sort_SE_library, \
                          args = (settings,library,mode,first_task)) \
                          for library in sorted(settings["libraries"])]
        pool.close()
        pool.join()
        for r in results:
            r.get()
            
    def sam_sort_PE_library(self,settings,library,mode,first_task):
        '''
        For paired-end data.
        Sorting one library at the time.
        '''
        print("\t"+library)
        #set file names
        if mode != "all":
            input_file = os.path.join(settings["--output"],"align","round_"+mode,\
                                      library+"_"+mode+".sam")
            sorted_file = os.path.join(settings["--output"],"sam_sort","round_"+mode,\
                                       library+"_"+mode+"_sort.sam")
            unmapped_file = os.path.join(settings["--output"],"sam_sort","round_"+mode,\
                                         library+"_"+mode+"_unmapped.sam")
        else:
            if first_task == "sam_sort":
                input_file = os.path.join(settings["--input"],library+".sam")
            else:
                input_file = os.path.join(settings["--output"],"align",\
                                      library+".sam")
            sorted_file = os.path.join(settings["--output"],"sam_sort",\
                                       library+"_sort.sam")
            unmapped_file = os.path.join(settings["--output"],"sam_sort",\
                                         library+"_unmapped.sam")

        #open files
        f_input = open(input_file)
        f_sorted = open(sorted_file,"w")
        f_unmapped = open(unmapped_file,"w")
        
        #parse input_file
        mate1_new = f_input.readline()
        ##pars and write header
        while mate1_new[0] == "@":
            if mode != "II": 
                f_sorted.write(mate1_new)
            f_unmapped.write(mate1_new)
            mate1_new = f_input.readline()
            if mate1_new == "":
                break
            
        ##statistics    
        mapped_reads = 0
        mapping_distribution = {}
        unmapped_reads = 0
        problem_reads = 0
        ## pars rest of the input file
        for line in f_input:
            #eliminating unmapped reads and read where mate is on mapped
            #(discordantly mapped reads)
            while mate1_new.split("\t")[2] =="*" or int(mate1_new.split("\t")[1]) & 8: 
                unmapped_reads += 1
                f_unmapped.write(mate1_new)
                f_unmapped.write(line)
                mate1_new = f_input.readline()
                line = f_input.readline()
                if len(mate1_new) == 0: #last lane
                    break
            #mapped reads    
            else:
                mappings = []
                mapping_number = 0
                read_id = mate1_new.split("\t")[0]
                mate1 = mate1_new
                mate2 = line
                #saving sequencing quality and orienting it to forward strain
                if int(mate1.split("\t")[1]) & 16: 
                    seq_qual1 = mate1.split("\t")[10][::-1]
                    seq_qual2 = mate2.split("\t")[10][::-1]
                else:
                    seq_qual1 = mate1.split("\t")[10]
                    seq_qual2 = mate2.split("\t")[10]
                    
                mapping_quality = int(mate1.split("\t")[11].split(":")[-1].strip())+\
                                  int(mate2.split("\t")[11].split(":")[-1].strip())
                mappings.append([mate1,mate2,mapping_quality]) # adding first pair to the list
                mate1_new = f_input.readline()
                # adding next pairs to the list, but eliminating discordantly mapped reads:
                while mate1_new.split("\t")[0] == read_id and \
                      not int(mate1.split("\t")[1]) & 8 and \
                      not int(mate1.split("\t")[1]) & 4:
                    mate2_new = f_input.readline()
                    mapping_quality = int(mate1_new.split("\t")[11].split(":")[-1].strip())+\
                                      int(mate2_new.split("\t")[11].split(":")[-1].strip())
                    mappings.append([mate1_new,mate2_new,mapping_quality]) 
                    mate1_new = f_input.readline()
                #looks up the best mapping quality 
                max_qual= -1000
                for pair in mappings:            
                    if max_qual < pair[2]:
                        max_qual = pair[2]
                        
                #looks up the mapping in best quality class
                better_mappings = []
                for pair in mappings: 
                    if max_qual == pair[2]:
                        better_mappings.append(pair[:-1])

                #removes overlapping mappings  
                best_mappings =[]
                wide_mappings = set()
                for pair in better_mappings:
                    #creates a set of positions where reads from pair were mapped
                    if int(pair[0].split("\t")[1]) & 16: 
                        set1 = set(range(int(pair[0].split("\t")[7]),\
                                         int(pair[0].split("\t")[3])+\
                                         abs(int(pair[0].split("\t")[5][:-1]))))
                    else:
                        set1 = set(range(int(pair[0].split("\t")[3]),\
                                         int(pair[0].split("\t")[7])+\
                                         abs(int(pair[1].split("\t")[5][:-1]))))
                        
                    #compare mappings among better_mappings group    
                    for pair2 in better_mappings: 
                        if pair != pair2: #do not compare thing with itself
                            #creates set of positions where reads from pair2 were mapped
                            if int(pair[0].split("\t")[1]) & 16: 
                                set2 = set(range(int(pair2[0].split("\t")[7]),\
                                                 int(pair2[0].split("\t")[3])+\
                                                 abs(int(pair2[0].split("\t")[5][:-1]))))
                            else:
                                set2 = set(range(int(pair2[0].split("\t")[3]),\
                                                 int(pair2[0].split("\t")[7])+\
                                                 abs(int(pair2[1].split("\t")[5][:-1]))))
                                
                            #this eliminates also weird pairs if some of the adapter is left
                            #(due to the wrong base calling) and one mate is mapping nearby
                            #region causing 3 mappings (one long + two short but overlaping)
                            if len(set1.intersection(set2)) >= \
                                        min(int(pair[0].split("\t")[5][:-1]),\
                                        int(pair2[0].split("\t")[5][:-1])):
                                #compares length of the mapping and longer goes to wide
                                if abs(int(pair2[0].split("\t")[8])) > \
                                   abs(int(pair[0].split("\t")[8])):
                                    wide_mappings.add(tuple(pair2))
                                else:
                                    wide_mappings.add(tuple(pair))
                                    break
                                
                #passing on only short mappings
                ##checking if there are wide mappings and
                ##allowing only those which are shorter
                if len(wide_mappings) > 0:  
                    for pair in better_mappings:
                        if tuple(pair) not in wide_mappings:
                            best_mappings.append(pair)
                            mapping_number += 1
                else:
                    best_mappings = better_mappings
                    mapping_number = len(best_mappings)

                #this situation happens when some adapter is left
                #and 3 mappings are generated (one long and 2 short).
                if len(best_mappings) == 0: 
                    problem_reads += 1
                    continue


                #editing seq quality field and samflag field
                ##first mapping
                read1 = best_mappings[0][0].split("\t")
                read2 = best_mappings[0][1].split("\t")
                ###seq quality is oriented in correct direction
                ###in case the first mapping was removed
                if int(read1[1]) & 16: #orienting quality to reverse strain
                    read1[10] = seq_qual1[::-1]
                    read2[10] = seq_qual2[::-1]
                else:
                    read1[10] = seq_qual1
                    read2[10] = seq_qual2
                ###secondary alignment is removed from sam flag
                ###in case the first mapping was removed
                read1[1] = int(read1[1])
                read2[1] = int(read2[1])
                if read1[1] & 256:
                    read1[1] = read1[1]-256
                if read2[1] & 256:
                    read2[1] = read2[1]-256
                read1[1] = str(read1[1])
                read2[1] = str(read2[1])        
                best_mappings[0][0] = "\t".join(read1)
                best_mappings[0][1] = "\t".join(read2)

                ###next mappings
                for x in range(1,len(best_mappings)):
                    read1 = best_mappings[x][0].split("\t")
                    read2 = best_mappings[x][1].split("\t")
                    ###seq quality is replaced for secondary alignment
                    ###(just in case)
                    read1[10] = "*"
                    read2[10] = "*"
                    ###secondary alignment is added to sam flag
                    ###if it was not there (just in case)
                    read1[1] = int(read1[1])
                    read2[1] = int(read2[1])
                    if read1[1] & 256:
                        pass
                    else:
                        read1[1] = read1[1]+256
                    if read2[1] & 256:
                        pass
                    else:
                        read2[1] = read2[1]+256
                    read1[1] = str(read1[1])
                    read2[1] = str(read2[1])
                    best_mappings[x][0] = "\t".join(read1)
                    best_mappings[x][1] = "\t".join(read2)

                #write selected mappings to the file                   
                for pair in best_mappings:
                    f_sorted.write(pair[0])
                    f_sorted.write(pair[1])

                #get mapping statistics
                mapped_reads += 1
                if mapping_number in mapping_distribution:
                    mapping_distribution[mapping_number] += 1
                else:
                    mapping_distribution[mapping_number] = 1

        #write mapping statistics
        self.write_sort_stat(settings,mode,library,mapped_reads,unmapped_reads,
                        problem_reads,mapping_distribution)
   
        f_input.close()
        f_sorted.close()
        f_unmapped.close()

    def sam_sort_SE_library(self,settings,library,mode,first_task):
        '''
        For single-end data
        Sorting one library at the time.
        '''
        print("\t"+library)
        #set file names
        if mode != "all":
            input_file = os.path.join(settings["--output"],"align","round_"+mode,\
                                      library+"_"+mode+".sam")
            sorted_file = os.path.join(settings["--output"],"sam_sort","round_"+mode,\
                                       library+"_"+mode+"_sort.sam")
            unmapped_file = os.path.join(settings["--output"],"sam_sort","round_"+mode,\
                                         library+"_"+mode+"_unmapped.sam")
        else:
            if first_task == "sam_sort":
                input_file = os.path.join(settings["--input"],library+".sam")
            else:
                input_file = os.path.join(settings["--output"],"align",\
                                      library+".sam")
            sorted_file = os.path.join(settings["--output"],"sam_sort",\
                                       library+"_sort.sam")
            unmapped_file = os.path.join(settings["--output"],"sam_sort",\
                                         library+"_unmapped.sam")
            
        #open files
        f_input = open(input_file)
        f_sorted = open(sorted_file,"w")
        f_unmapped = open(unmapped_file,"w")
        
        #pars input_file
        mate1_new = f_input.readline()
        ##pars and write header
        while mate1_new[0] == "@":
            if mode != "II": 
                f_sorted.write(mate1_new)
            f_unmapped.write(mate1_new)
            mate1_new = f_input.readline()
            if mate1_new == "":
                break
            
        ##statistics    
        mapped_reads = 0
        mapping_distribution = {}
        unmapped_reads = 0
        problem_reads = 0

        ## pars rest of the input file
        for line in f_input:
            #eliminating unmapped reads
            while mate1_new.split("\t")[2] =="*" and int(mate1_new.split("\t")[1]) & 4: 
                unmapped_reads += 1
                f_unmapped.write(mate1_new)
                mate1_new = line
                line = f_input.readline()
                if len(mate1_new) == 0: #last lane
                    break
            #mapped reads    
            else:
                mappings = []
                mapping_number = 0
                read_id = mate1_new.split("\t")[0]
                mate1 = mate1_new
                #saving sequencing quality and orienting it to forward strain
                if int(mate1.split("\t")[1]) & 16: 
                    seq_qual1 = mate1.split("\t")[10][::-1]
                else:
                    seq_qual1 = mate1.split("\t")[10]

                mapping_quality = int(mate1.split("\t")[11].split(":")[-1].strip())
                mappings.append([mate1,mapping_quality]) # adding first pair to the list
                mate1_new = line
                while mate1_new.split("\t")[0] == read_id: # adding next pairs to the list
                    mapping_quality = int(mate1_new.split("\t")[11].split(":")[-1].strip())
                    mappings.append([mate1_new,mapping_quality])
                    mate1_new = f_input.readline()
                #looks up the best mapping quality 
                max_qual= -1000
                for pair in mappings:            
                    if max_qual < pair[1]:
                        max_qual = pair[1]
                        
                #looks up the mapping in best quality class
                better_mappings = []
                for pair in mappings: 
                    if max_qual == pair[1]:
                        better_mappings.append(pair[:-1])



                #editing seq quality field and samflag field
                ##first mapping
                read1 = better_mappings[0][0].split("\t")

                ###seq quality is oriented in correct direction
                ###in case the first mapping was removed
                if int(read1[1]) & 16: #orienting quality to reverse strain
                    read1[10] = seq_qual1[::-1]
                else:
                    read1[10] = seq_qual1

                ###secondary alignment is removed from sam flag
                ###in case the first mapping was removed
                read1[1] = int(read1[1])
                if read1[1] & 256:
                    read1[1] = read1[1]-256
                read1[1] = str(read1[1])       
                better_mappings[0][0] = "\t".join(read1)

                ###next mappings
                for x in range(1,len(better_mappings)):
                    read1 = better_mappings[x][0].split("\t")

                    ###seq quality is replaced for secondary alignment
                    ###(just in case)
                    read1[10] = "*"

                    ###secondary alignment is added to sam flag
                    ###if it was not there (just in case)
                    read1[1] = int(read1[1])
                    if read1[1] & 256:
                        pass
                    else:
                        read1[1] = read1[1]+256
                    read1[1] = str(read1[1])
                    better_mappings[x][0] = "\t".join(read1)

                #write selected mappings to the file                   
                for pair in better_mappings:
                    f_sorted.write(pair[0])

                #get mapping statistics
                mapped_reads += 1
                mapping_number = len(better_mappings)
                if mapping_number in mapping_distribution:
                    mapping_distribution[mapping_number] += 1
                else:
                    mapping_distribution[mapping_number] = 1

        #write mapping statistics
        self.write_sort_stat(settings,mode,library,mapped_reads,unmapped_reads,
                        problem_reads,mapping_distribution)
            
        f_input.close()
        f_sorted.close()
        f_unmapped.close()

    def write_sort_stat(self,settings,mode,library,mapped_reads,unmapped_reads,
                        problem_reads,mapping_distribution):
        '''
        Write sorting data to file.
        '''
        if mode != "all":
            report_file = os.path.join(settings["--output"],"sam_sort","round_"+mode,\
                                       "sam_sort_info",library+"_"+mode+"sam_sortinfo.log")

        else:
            report_file = os.path.join(settings["--output"],"sam_sort",\
                                       "sam_sort_info",library+"_sam_sortinfo.log")
        f_report = open(report_file,"w")
        f_report.write("Mapped reads "+str(mapped_reads) +"\n")
        f_report.write("Unmapped reads "+str(unmapped_reads) +"\n")
        f_report.write("Problematic (with adapter) overlapping reads "+str(problem_reads) +"\n")
        f_report.write("Mapping distribution:" + "\n")
        for element in sorted(mapping_distribution):
            f_report.write(str(element) +"\t" + str(mapping_distribution [element])+ "\n")
        f_report.close()


    def sam_unmapped_to_fq_PE(self,settings):
        '''
        FOr paired-end reads.
        Converts sam file to fastq file.
        '''
        for library in sorted(settings["libraries"]):           
            #open files
            f_in = open(os.path.join(settings["--output"],"sam_sort","round_I",\
                                     library+"_I_unmapped.sam"))
            f_out1 = open(os.path.join(settings["--output"],"align","round_I",\
                                     library+"_unmapped_1.fq"), "w")
            f_out2 = open(os.path.join(settings["--output"],"align","round_I",\
                                     library+"_unmapped_2.fq"), "w")
            for line in f_in:
                #remove header
                while line.startswith("@"):
                    line = f_in.readline()
                #stop at file end
                if line == "":
                    continue
                seq_id = line.split("\t")[0]
                seq1 = line.split("\t")[9]
                seq_qual1=line.split("\t")[10]
                line = f_in.readline()
                seq2 = line.split("\t")[9]
                seq_qual2=line.split("\t")[10]
                f_out1.write("@"+seq_id+"/1\n")
                f_out1.write(seq1+"\n")
                f_out1.write("+\n")
                f_out1.write(seq_qual1+"\n")
                f_out2.write("@"+seq_id+"/2\n")
                f_out2.write(seq2+"\n")
                f_out2.write("+\n")
                f_out2.write(seq_qual2+"\n")

            f_in.close()
            f_out1.close()
            f_out2.close()

    def sam_unmapped_to_fq_SE(self,settings):
        '''
        For single-end reads.
        Converts sam file to fastq file.
        '''
        for library in sorted(settings["libraries"]):          
            #open files
            f_in = open(os.path.join(settings["--output"],"sam_sort","round_I",\
                                     library+"_I_unmapped.sam"))
            f_out1 = open(os.path.join(settings["--output"],"align","round_I",\
                                     library+"_unmapped.fq"), "w")
            for line in f_in:
                #remove header
                while line.startswith("@"):
                    line = f_in.readline()
                #stop at file end
                if line == "":
                    continue
                seq_id = line.split("\t")[0]
                seq1 = line.split("\t")[9]
                seq_qual1=line.split("\t")[10]
                f_out1.write("@"+seq_id+"/1\n")
                f_out1.write(seq1+"\n")
                f_out1.write("+\n")
                f_out1.write(seq_qual1+"\n")

            f_in.close()
            f_out1.close()

    def join_sam(self,settings):
        '''
        Joins SAM files from round I and II.
        In addition creates general info file
        '''
        for library in sorted(settings["libraries"]):
            #join sam files
            join_sam_command = (
                "cat",
                os.path.join(settings["--output"],"sam_sort","round_I",library+"_I_sort.sam"),
                os.path.join(settings["--output"],"sam_sort","round_II",library+"_II_sort.sam"),
                ">", os.path.join(settings["--output"],"sam_sort",library+"_sort.sam")
                )
            os.system(" ".join(list(join_sam_command)))

            #create common info file
            report_file1 = os.path.join(settings["--output"],"sam_sort","round_I",\
                                           "sam_sort_info",library+"_I_sam_sortinfo.log")
            report_file2 = os.path.join(settings["--output"],"sam_sort","round_II",\
                                           "sam_sort_info",library+"_II_sam_sortinfo.log")
            ##read in from file I
            with open(report_file1) as f_in1:
                mapping_distribution = {}
                for line in f_in1:
                    if line.startswith("Mapped"):
                        mapped_reads = int(line.strip().split(" ")[2])
                    if line.startswith("Problematic"):
                        problem_reads = int(line.strip().split(" ")[3])
                    elif len(line.strip().split("\t")) == 2:
                        mapping_distribution[line.strip().split("\t")[0]] = \
                                                            int(line.strip().split("\t")[1])
            ##read in from file II
            with open(report_file2) as f_in2:
                for line in f_in2:
                    if line.startswith("Mapped"):
                        mapped_reads += int(line.strip().split(" ")[2])
                    elif line.startswith("Unmapped"):
                        unmapped_reads = int(line.strip().split(" ")[2])
                    elif line.startswith("Problematic"):
                        problem_reads += int(line.strip().split(" ")[3])
                    elif len(line.strip().split("\t")) == 2:
                        if line.strip().split("\t")[0] not in mapping_distribution:
                            mapping_distribution[line.strip().split("\t")[0]] = \
                                                            int(line.strip().split("\t")[1])
                        else:
                            mapping_distribution[line.strip().split("\t")[0]] += \
                                                            int(line.strip().split("\t")[1])
            ##write new reportfile
            report_file = os.path.join(settings["--output"],"sam_sort",\
                                           "sam_sort_info",library+"_aligninfo.log")                
            with open(report_file, "w") as f_out:
                f_out.write("Mapped reads "+str(mapped_reads) +"\n")
                f_out.write("Unmapped reads "+str(unmapped_reads) +"\n")
                if settings["paired"]:
                    f_out.write("Problematic overlapping reads "+str(problem_reads) +"\n")
                f_out.write("Mapping distribution" + "\n")
                for element in sorted(mapping_distribution):
                    f_out.write(str(element) +"\t" + str(mapping_distribution [element])+ "\n")
                                     
    def test_errors(self,settings,mode):
        '''
        Testing existence and errors in info files
        '''
        if mode != "all":
            info_folder = os.path.join(settings["--output"],"sam_sort",\
                                       "round_"+mode,"sam_sort_info")
        else:
            info_folder = os.path.join(settings["--output"],"sam_sort","sam_sort_info")
        #testing log file
        if settings["paired"]:
            wanted_read_beginnings = {"Mapped","Unmapped","Problematic","Mapping"}
        else:
            wanted_read_beginnings = {"Mapped","Unmapped","Mapping"}
        if len(os.listdir(info_folder)) == 0:
            sys.exit('Task "sort" incomplete! Infofiles missing in folder ' + info_folder)        
        for file in os.listdir(info_folder):
            with open(os.path.join(info_folder,file)) as f_in:
                read_beginnings = set()
                for line in f_in:
                    read_beginnings.add(line.split(" ")[0].split("\t")[0])
                if not wanted_read_beginnings.issubset(read_beginnings):
                    sys.exit('Error in task "sort", check infofile: '+\
                             os.path.join(info_folder,file))
