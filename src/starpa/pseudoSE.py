#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
1)Convert paired-end reads to pseudo single-end reads.
2)Eliminate reads/alignments with too many mismatches
3)Eliminate reads with too many alignments
4)Add NH:i: tag for count alignment number

To do:

'''

import os
import sys
import copy
from itertools import groupby
import multiprocessing as mp

from pyfaidx import Fasta

class pseudoSE():
    def __init__(self,settings,first_task):
        self.settings = settings
        self.make_folder(settings)
        if settings["paired"]:
            self.make_pseudoSE(settings,first_task)
        else:
            self.make_pseudoSE_SE(settings,first_task)
        self.test_errors(settings)

    def make_folder(self,settings):
        '''
        Make output folder for task if needed
        '''
        if not os.path.exists(os.path.join(settings["--output"],"pseudoSE")):
            os.makedirs(os.path.join(settings["--output"],"pseudoSE"))
        if not os.path.exists(os.path.join(settings["--output"],"pseudoSE","mismatched")):
            os.makedirs(os.path.join(settings["--output"],"pseudoSE","mismatched"))
        if not os.path.exists(os.path.join(settings["--output"],"pseudoSE","too_many_matches")):
            os.makedirs(os.path.join(settings["--output"],"pseudoSE","too_many_matches"))
        if not os.path.exists(os.path.join(settings["--output"],"pseudoSE","pseudoSE_info")):
            os.makedirs(os.path.join(settings["--output"],"pseudoSE","pseudoSE_info"))
        if settings["pseudoSE"]["pseudoSE_oligoA"]:
            if not os.path.exists(os.path.join(settings["--output"],"pseudoSE","oligoA")):
                os.makedirs(os.path.join(settings["--output"],"pseudoSE","oligoA"))
                
    def make_pseudoSE(self,settings,first_task):
        '''
        Generate pseudo single-end reads from paired-end reads
        '''
            
        pool = mp.Pool(processes=settings["CPUs"])
        results = [pool.apply_async(self.make_pseudoSE_library, \
                          args = (settings,library,first_task)) \
                          for library in sorted(settings["libraries"])]
        pool.close()
        pool.join()
        for r in results:
            r.get()

    def make_pseudoSE_SE(self,settings,first_task):
        '''
        Filter single-end reads suitable from furtehr tasks
        '''
        pool = mp.Pool(processes=settings["CPUs"])
        results = [pool.apply_async(self.make_pseudoSE_SE_library, \
                          args = (settings,library,first_task)) \
                          for library in sorted(settings["libraries"])]
        pool.close()
        pool.join()
        for r in results:
            r.get()

    def make_pseudoSE_library(self,settings,library,first_task):
        '''
        Create pseudo single-end reads from paired-end reads of one library
        '''
        print("\t"+library)
        genome = Fasta(settings["genome"],one_based_attributes=False)
        #set filenames
        if first_task == "pseudoSE":
            input_file = os.path.join(settings["--input"],library+\
                                        settings["pseudoSE"]["pseudoSE_input_file_suffix"])
        else:
            input_file = os.path.join(settings["--output"],"sam_sort",\
                                          library+"_sort.sam")
        pseudoSE_file = os.path.join(settings["--output"],"pseudoSE",\
                                          library+"_pseudoSE.sam")
        #reads with too many mismatches:
        mismatch_file = os.path.join(settings["--output"],"pseudoSE","mismatched",\
                                          library+"_pseudoSE_mismatch.sam")
        #reads with too many mappings:
        many_match_file = os.path.join(settings["--output"],"pseudoSE","too_many_matches",\
                                          library+"_pseudoSE_multimatch.sam")

        if settings["pseudoSE"]["pseudoSE_oligoA"]:
            #mappings which have oligoA in the end:
            oligoA_file = os.path.join(settings["--output"],"pseudoSE","oligoA",\
                                              library+"-oligoA_pseudoSE.sam")
            #mappings which whould be discarded with oligoA,
            #but not without oligoA:
            oligoA_mismatch_file = os.path.join(settings["--output"],"pseudoSE","oligoA",\
                                              library+"-oligoA-mm_pseudoSE.sam")

        #open files
        f_input = open(input_file)
        f_pseudoSE = open(pseudoSE_file,"w")
        f_mismatch = open(mismatch_file,"w")
        f_many_match = open(many_match_file,"w")
        if settings["pseudoSE"]["pseudoSE_oligoA"]:
            f_oligoA = open(oligoA_file,"w") #all oligoA matches
            #oligoA matches which would not have passed when oligoA would not be allowed:
            f_oligoA_mismatch = open(oligoA_mismatch_file,"w") 

        #statistics
        total_reads = 0
        total_mappings = 0
        mapped_reads = 0
        mappings_counter = 0
        many_match_reads = 0
        many_match_mappings = 0
        mismatched_reads = 0
        mismatched_mappings = 0
        #will be used only if oligoA is considered
        oligoA_reads = 0 #oligoA reads where all mappings are oligoA
        oligoA_mappings = 0 #all oligoA mappings
        #oligoA reads where all mappings are oligoA mismatches:
        oligoA_mismatch_reads = 0 
        #all oligoA mappings which would have not passed when oligoA
        #would not be allowed:
        oligoA_mismatch_mappings = 0
        
        mapping_distribution = {}

        #pars input_file
        mate1_new = f_input.readline()
        ##pars and write header
        while mate1_new[0] == "@":
            f_pseudoSE.write(mate1_new)
            f_mismatch.write(mate1_new)
            f_many_match.write(mate1_new)
            if settings["pseudoSE"]["pseudoSE_oligoA"]:
                f_oligoA.write(mate1_new)
                f_oligoA_mismatch.write(mate1_new)
            mate1_new = f_input.readline()
            if mate1_new == "": #end of the file
                break

        ## pars rest of the input file
        for line in f_input:
                
            #skip unmapped reads (just in case)
            if int(mate1_new.split("\t")[1]) & 4: #unmapped read 
                mate1_new = f_input.readline()
                continue
            
            #if mate unmapped skip pair (just in case)
            if int(line.split("\t")[1]) & 4: #unmapped mate
                mate1_new = f_input.readline()                    
                continue
            
            #getting all mappings of single read
            ##getting first pair
            mappings = []
            read_id = mate1_new.strip().split("\t")[0]
            mate1 = mate1_new.strip().split("\t")
            mate2 = line.strip().split("\t")
            mappings.append([mate1,mate2]) 
            ##adding next pairs to the list
            mate1_new = f_input.readline()
            while mate1_new.split("\t")[0] == read_id:
                mate2_new = f_input.readline()
                mappings.append([mate1_new.strip().split("\t"),mate2_new.strip().split("\t")])
                mate1_new = f_input.readline()

            total_reads += 1
            total_mappings += len(mappings)

            #Convert PE to SE
            ##set raw quality from mapping which has quality by scanning all mappings
            for mapping in mappings:
                if mapping[0][10] != "*":
                    if int(mappings[0][0][1]) & 16: # orienting quality to forward strand
                        raw_qual1 = mappings[0][0][10][::-1]
                        raw_qual2 = mappings[0][1][10][::-1]
                    else:
                        raw_qual1 = mappings[0][0][10]
                        raw_qual2 = mappings[0][1][10]
                    break
            read_lengths = []
            for mapping in mappings:

                ##get quality mapping
                ##secondary mappings which do not have seq quality are provided quality
                ##from primary mapping
                if mapping[0][10] == "*": #
                    if int(mapping[0][1]) & 16: # orienting quality to reverse strand
                        qual1 = raw_qual1[::-1]
                        qual2 = raw_qual2[::-1]
                    else:
                        qual1 = raw_qual1
                        qual2 = raw_qual2                    
                else:
##                    if len(mapping) == 1:
##                        print(mapping)
                    #print(mapping)
                    if int(mapping[0][1]) & 16: # orienting quality to reverse strand
                        qual1 = mapping[0][10][::-1]
                        qual2 = mapping[1][10][::-1]
                    else:
                        qual1 = mapping[0][10]
                        qual2 = mapping[1][10]                         
                ##if needed convert seq quality to phred+33
                if settings["pseudoSE"]["pseudoSE_quality_base"] == 64:
                    qual1 = self.convert_to_phred33(qual1)
                    qual2 = self.convert_to_phred33(qual2)

                ##combine sequence and quality
                seq1 =  mapping[0][9]
                seq2 =  mapping[1][9]
                read_lengths.append([len(seq1),len(seq2)])
                seq1pos1 = int(mapping[0][3])
                seq1pos2 = int(mapping[0][3]) + int(mapping[0][5][:-1])
                seq2pos1 = int(mapping[1][3])
                seq2pos2 = int(mapping[1][3]) + int(mapping[1][5][:-1])
                chrom = str(mapping[0][2])
                if seq1pos1 > seq2pos2: #no overlap
                    #full sequence:
                    mapping[0][9] = seq2+genome[chrom][seq2pos2-1:seq1pos1-1].seq+seq1 
                    mapping[0][10] = qual2+"J"*(seq1pos1-seq2pos2)+qual1
                elif seq2pos1 > seq1pos2: #no overlap
                    #full sequence:
                    mapping[0][9] = seq1+genome[chrom][seq1pos2-1:seq2pos1-1].seq+seq2
                    mapping[0][10] = qual1+"J"*(seq2pos1-seq1pos2)+qual2
                else: #some or full overlap
                    offset1 = seq1pos1-seq2pos1
                    offset2 = seq1pos2-seq2pos2
                    #print(seq1,seq2,qual1,qual2,offset1,offset2)
                    mapping[0][9],mapping[0][10] = \
                                    self.best_seq(seq1,seq2,qual1,qual2,offset1,offset2)

                ##sam fields from PE to SE
                mapping[0][3] = str(min(int(mapping[0][3]),int(mapping[0][7])))
                mapping[0][6],mapping[0][7],mapping[0][8] = "*",str(0),str(0)    
                #create new CIGAR string
                mapping[0][5] = str(len(mapping[0][9]))+"M"
                
            #Check mismatches
            good_mappings = []
            oligoA_passed_mismatches = 0
            for z,mapping in enumerate(mappings):
                seq = mapping[0][9]
                ref = genome[chrom][int(mapping[0][3])-1:\
                                    int(mapping[0][3])+len(mapping[0][9])-1].seq
                #creating alignement, and counting matches and mismatches
                alignment,tot_mismatch = self.counting_mismatches(seq,ref)

                #adding tag with mismatching position
                mapping[0] = self.change_add_extra_tag("MD:Z:",alignment,mapping[0])

                #setting allowed mismatched to be depenent from the length of the sequence
                #and given number
                mismatch_limit = max(settings["pseudoSE"]["pseudoSE_allowed_mismatch"],\
                        ##take minimum from the two:
                        ###set precentage for the maximum sequenced are
                        ##min(settings["pseudoSE"]["pseudoSE_max_read_length"]*2*\
                        min(sum(read_lengths[z])*\
                            settings["pseudoSE"]["pseudoSE_mismatch_precentage"]/100,\
                            ###precentage from the full length SE sequence
                            len(seq)*settings["pseudoSE"]["pseudoSE_mismatch_precentage"]/100))
                
                #Testing mismatch levels   
                if tot_mismatch > mismatch_limit:
                    i = 1
                    if int(mapping[0][1]) & 16: #mapped to the negative strand
                        #if oligoA reads are included
                        if settings["pseudoSE"]["pseudoSE_oligoA"]:
                            #testing mismatches without oligoA
                            #and mapping is preserved if mismatches are below threshold
                            while seq[i-1] in {"T","N"}:
                                alignement,tot_mismatch = \
                                                self.counting_mismatches(seq[i+1:],ref[i+1:])
                                i +=1
                                #read considered as not haveing too many mismatches
                                if tot_mismatch <= mismatch_limit:
                                    #changing flag just to indicate strand and primary/secondary
                                    mapping[0] = self.change_sam_flag(mapping[0],good_mappings)
                                    f_oligoA_mismatch.write("\t".join(mapping[0])+"\n")
                                    oligoA_mismatch_mappings += 1
                                    oligoA_passed_mismatches += 1
                                    good_mappings.append(mapping[0])
                                    break
                            else: #too many mismatches
                                f_mismatch.write("\t".join(mapping[0])+"\n")
                                mismatched_mappings += 1
                             
                        else: #read considered with too many mismatches
                            f_mismatch.write("\t".join(mapping[0])+"\n")
                            mismatched_mappings += 1

                    else: # mapped on positive strand
                         #if oligoA reads are included
                        if settings["pseudoSE"]["pseudoSE_oligoA"]:
                            #testing mismatches without oligoA
                            #and mapping is preserved if mismatches are below threshold
                            while seq[-i] in {"A","N"}:
                                alignement,tot_mismatch = \
                                                    self.counting_mismatches(seq[:-i],ref[:-i])
                                i +=1
                                #read considered as not haveing too many mismatches
                                if tot_mismatch <= mismatch_limit:
                                    #changing flag just to indicate strand and primary/secondary
                                    mapping[0] = self.change_sam_flag(mapping[0],good_mappings)
                                    f_oligoA_mismatch.write("\t".join(mapping[0])+"\n")
                                    oligoA_mismatch_mappings += 1
                                    oligoA_passed_mismatches += 1
                                    good_mappings.append(mapping[0])
                                    break
                            else: #too many mismatches
                                f_mismatch.write("\t".join(mapping[0])+"\n")
                                mismatched_mappings += 1
                             
                        else: #read considered with too many mismatches
                            f_mismatch.write("\t".join(mapping[0])+"\n")
                            mismatched_mappings += 1

                else: #mismatches below threshold
                    #changing flag just to indicate strand and primary/secondary
                    mapping[0] = self.change_sam_flag(mapping[0],good_mappings)
                    good_mappings.append(mapping[0])

            
            #if all reads are mismatched
            if len(good_mappings) == 0:
                mismatched_reads += 1
                continue
            
            #cound oligoA mismatch reads (reads which would fail if oligoA is not considered
            if oligoA_passed_mismatches == len(good_mappings):
                oligoA_mismatch_reads += 1
                
            #Check number of mappings
            ##add tag on match number
            for mapping in good_mappings:
                mapping = self.change_add_extra_tag("NH:i:",len(good_mappings),mapping)

            ##skip read with many mapping and write rest to the file
            if len(good_mappings) > settings["pseudoSE"]["pseudoSE_max_mappings"]:
                for mapping in good_mappings:
                    f_many_match.write("\t".join(mapping)+"\n")
                many_match_reads += 1
                many_match_mappings += len(good_mappings)
                continue
            else:
                for mapping in good_mappings:
                    f_pseudoSE.write("\t".join(mapping)+"\n")
                mapped_reads += 1
                #if len(good_mappings) > 1:
                    #print(good_mappings)
                mappings_counter += len(good_mappings)                    
                if len(good_mappings) not in mapping_distribution:
                    mapping_distribution[len(good_mappings)] = 1
                else:
                    mapping_distribution[len(good_mappings)] += 1

                #saving reads separately which have mismatch A at 3'
                if settings["pseudoSE"]["pseudoSE_oligoA"]:
                    oligoA = 0
                    for mapping in good_mappings:
                        ref = genome[mapping[2]][int(mapping[3])-1:\
                                                 int(mapping[3])-1+\
                                                 int(mapping[5].strip("M"))].seq
                        #all reads which do have A-s or N-s in the 3' end and
                        #a single mismatch in this region are considered as
                        #reads with oligoA tail
                        i = 1
                        if int(mapping[1]) & 16:
                            while mapping[9][i-1] in {"T","N"}:
                                if ref[i-1] != "T":
                                    f_oligoA.write("\t".join(mapping)+"\n")
                                    oligoA_mappings += 1
                                    oligoA +=1
                                    break
                                i +=1
                        else:
                            while mapping[9][-i] in {"A","N"}:
                                if ref[-i] != "A":
                                    f_oligoA.write("\t".join(mapping)+"\n")
                                    oligoA_mappings += 1
                                    oligoA +=1
                                    break
                                i +=1
                    #count oligoA reads
                    if oligoA == len(good_mappings):
                        oligoA_reads += 1
                            
        #write info file
        self.write_statistics(settings,library,total_reads,mismatched_reads,many_match_reads,
                              mapped_reads,oligoA_reads,oligoA_mismatch_reads,total_mappings,
                              mismatched_mappings,many_match_mappings,mappings_counter,
                              oligoA_mappings,oligoA_mismatch_mappings,mapping_distribution)
        #close files
        f_input.close()
        f_pseudoSE.close()
        f_mismatch.close()
        f_many_match.close()
        if settings["pseudoSE"]["pseudoSE_oligoA"]:
            f_oligoA.close()
            f_oligoA_mismatch.close()

    def make_pseudoSE_SE_library(self,settings,library,first_task):
        '''
        Filter single-end reads suitable from furtehr tasks from one library.
        '''
        genome = Fasta(settings["genome"],one_based_attributes=False)
        print("\t"+library)
        #set filenames
        if first_task == "pseudoSE":
            input_file = os.path.join(settings["--input"],library+\
                                        settings["pseudoSE"]["pseudoSE_input_file_suffix"])
        else:
            input_file = os.path.join(settings["--output"],"sam_sort",\
                                          library+"_sort.sam")
        pseudoSE_file = os.path.join(settings["--output"],"pseudoSE",\
                                          library+"_pseudoSE.sam")
        #reads with too many mismatches:
        mismatch_file = os.path.join(settings["--output"],"pseudoSE","mismatched",\
                                          library+"_pseudoSE_mismatch.sam")
        #reads with too many mappings:
        many_match_file = os.path.join(settings["--output"],"pseudoSE","too_many_matches",\
                                          library+"_pseudoSE_multimatch.sam")
        info_file = os.path.join(settings["--output"],"pseudoSE","pseudoSE_info",\
                                          library+"_pseudoSEinfo.log")
        if settings["pseudoSE"]["pseudoSE_oligoA"]:
            #mappings which have oligoA in the end:
            oligoA_file = os.path.join(settings["--output"],"pseudoSE","oligoA",\
                                              library+"-oligoA_pseudoSE.sam")
            #mappings which whould be discarded with oligoA,
            #but not without oligoA:
            oligoA_mismatch_file = os.path.join(settings["--output"],"pseudoSE","oligoA",\
                                              library+"-oligoA-mm_pseudoSE.sam")
        else:
            oligoA_file = ""
            oligoA_mismatch_file = ""

        #open files
        f_input = open(input_file)
        f_pseudoSE = open(pseudoSE_file,"w")
        f_mismatch = open(mismatch_file,"w")
        f_many_match = open(many_match_file,"w")
        f_info = open(info_file,"w")
        if settings["pseudoSE"]["pseudoSE_oligoA"]:
            f_oligoA = open(oligoA_file,"w") #all oligoA matches
            #oligoA matches which would not have passed when oligoA would not be allowed:
            f_oligoA_mismatch = open(oligoA_mismatch_file,"w")
        else:
            f_oligoA = ""
            f_oligoA_mismatch = ""

        #statistics
        total_reads = 0
        total_mappings = 0
        mapped_reads = 0
        mappings_counter = 0
        many_match_reads = 0
        many_match_mappings = 0
        mismatched_reads = 0
        mismatched_mappings = 0
        #is used only if oligA is considered
        oligoA_reads = 0 #oligoA reads where all mappings are oligoA
        oligoA_mappings = 0 #all oligoA mappings
        #oligoA reads where all mappings are oligoA mismatches:
        oligoA_mismatch_reads = 0 
        #all oligoA mappings which would have not passed when oligoA
        #would not be allowed:
        oligoA_mismatch_mappings = 0
            
        mapping_distribution = {}

        ##parse and write header
        mate1_new = f_input.readline()
        while mate1_new[0] == "@":
            f_pseudoSE.write(mate1_new)
            f_mismatch.write(mate1_new)
            f_many_match.write(mate1_new)
            if settings["pseudoSE"]["pseudoSE_oligoA"]:
                f_oligoA.write(mate1_new)
                f_oligoA_mismatch.write(mate1_new)
            mate1_new = f_input.readline()
            if mate1_new == "": #end of the file
                break


        ## pars rest of the input file
        for mate1_next in f_input:
    
            #skip unmapped reads (just in case)
            if int(mate1_new.split("\t")[1]) & 4: #unmapped read
                mate1_new = copy.deepcopy(mate1_next)
                continue
            
            #getting all mappings of single read
            ##getting first pair
            mappings = []
            read_id = mate1_new.strip().split("\t")[0]
            mate1 = mate1_new.strip().split("\t")
#                mate2 = line.strip().split("\t")
            mappings.append([mate1])
#                mappings.append([mate1,mate2]) 
            ##adding next pairs to the list
            while mate1_next.split("\t")[0] == read_id:
#                    mate2_new = f_input.readline()
#                    mappings.append([mate1_new.strip().split("\t"),mate2_new.strip().split("\t")])
                mappings.append([mate1_next.strip().split("\t")])
                mate1_next = f_input.readline()
                if mate1_next == "":
                    break
            mate1_new = copy.deepcopy(mate1_next)

            total_reads,total_mappings,mapped_reads,mappings_counter,\
            many_match_reads,many_match_mappings,mismatched_reads,mismatched_mappings,\
            oligoA_reads,oligoA_mappings,oligoA_mismatch_reads,oligoA_mismatch_mappings,\
            mapping_distribution = self.process_mappings_SE(settings,mappings,genome,
                                        f_input,f_pseudoSE,f_mismatch,f_many_match,f_info,f_oligoA,
                                        f_oligoA_mismatch,total_reads,total_mappings,mapped_reads,mappings_counter,
                                        many_match_reads,many_match_mappings,mismatched_reads,mismatched_mappings,
                                        oligoA_reads,oligoA_mappings,oligoA_mismatch_reads,oligoA_mismatch_mappings,
                                        mapping_distribution)
            print("C",mapping_distribution)
        #process last single mapping
        if mate1_new != "":
            mappings = [[mate1_next.strip().split("\t")]]
            total_reads,total_mappings,mapped_reads,mappings_counter,\
            many_match_reads,many_match_mappings,mismatched_reads,mismatched_mappings,\
            oligoA_reads,oligoA_mappings,oligoA_mismatch_reads,oligoA_mismatch_mappings,\
            mapping_distribution = self.process_mappings_SE(settings,mappings,genome,
                                        f_input,f_pseudoSE,f_mismatch,f_many_match,f_info,f_oligoA,
                                        f_oligoA_mismatch,total_reads,total_mappings,mapped_reads,mappings_counter,
                                        many_match_reads,many_match_mappings,mismatched_reads,mismatched_mappings,
                                        oligoA_reads,oligoA_mappings,oligoA_mismatch_reads,oligoA_mismatch_mappings,
                                        mapping_distribution)


        self.write_statistics(settings,library,total_reads,mismatched_reads,many_match_reads,
                              mapped_reads,oligoA_reads,oligoA_mismatch_reads,total_mappings,
                              mismatched_mappings,many_match_mappings,mappings_counter,
                              oligoA_mappings,oligoA_mismatch_mappings,mapping_distribution)                                

        f_input.close()
        f_pseudoSE.close()
        f_mismatch.close()
        f_many_match.close()
        f_info.close()
        if settings["pseudoSE"]["pseudoSE_oligoA"]:
            f_oligoA.close()
            #oligoA matches which would not have passed when oligoA would not be allowed:
            f_oligoA_mismatch.close() 


    def process_mappings_SE(self,settings,mappings,genome,
                            f_input,f_pseudoSE,f_mismatch,f_many_match,f_info,f_oligoA,
                            f_oligoA_mismatch,total_reads,total_mappings,mapped_reads,mappings_counter,
                            many_match_reads,many_match_mappings,mismatched_reads,mismatched_mappings,
                            oligoA_reads,oligoA_mappings,oligoA_mismatch_reads,oligoA_mismatch_mappings,
                            mapping_distribution):
        '''
        Processing mappings
        '''
        total_reads += 1
        total_mappings += len(mappings)

        #Convert PE to SE
        ##set raw quality from mapping which has quality by scanning all mappings
        for mapping in mappings:
            if mapping[0][10] != "*":
                if int(mappings[0][0][1]) & 16: # orienting quality to forward strand
                    raw_qual1 = mappings[0][0][10][::-1]
#                            raw_qual2 = mappings[0][1][10][::-1]
                else:
                    raw_qual1 = mappings[0][0][10]
#                            raw_qual2 = mappings[0][1][10]
                break
        for mapping in mappings:

            ##get quality mapping
            ##secondary mappings which do not have seq quality are provided quality
            ##from primary mapping
            if mapping[0][10] == "*": #
                if int(mapping[0][1]) & 16: # orienting quality to reverse strand
                    qual1 = raw_qual1[::-1]
#                            qual2 = raw_qual2[::-1]
                else:
                    qual1 = raw_qual1
#                            qual2 = raw_qual2                    
            else:
##                    if len(mapping) == 1:
##                        print(mapping)
                #print(mapping)
                if int(mapping[0][1]) & 16: # orienting quality to reverse strand
                    qual1 = mapping[0][10][::-1]
#                            qual2 = mapping[1][10][::-1]
                else:
                    qual1 = mapping[0][10]
#                            qual2 = mapping[1][10]
            ##if needed convert seq quality to phred+33
            if settings["pseudoSE"]["pseudoSE_quality_base"] == 64:
                qual1 = self.convert_to_phred33(qual1)
#                        qual2 = self.convert_to_phred33(qual2)

##                    ##combine sequence and quality
##                    seq1 =  mapping[0][9]
##                    seq2 =  mapping[1][9]
##                    seq1pos1 = int(mapping[0][3])
##                    seq1pos2 = int(mapping[0][3]) + int(mapping[0][5][:-1])
##                    seq2pos1 = int(mapping[1][3])
##                    seq2pos2 = int(mapping[1][3]) + int(mapping[1][5][:-1])
            chrom = str(mapping[0][2])
##                    if seq1pos1 > seq2pos2: #no overlap
##                        #full sequence:
##                        mapping[0][9] = seq2+genome[chrom][seq2pos2-1:seq1pos1-1]+seq1 
##                        mapping[0][10] = qual2+"J"*(seq1pos1-seq2pos2)+qual1
##                    elif seq2pos1 > seq1pos2: #no overlap
##                        #full sequence:
##                        mapping[0][9] = seq1+genome[chrom][seq1pos2-1:seq2pos1-1]+seq2
##                        mapping[0][10] = qual1+"J"*(seq2pos1-seq1pos2)+qual2
##                    else: #some or full overlap
##                        offset1 = seq1pos1-seq2pos1
##                        offset2 = seq1pos2-seq2pos2
##                        #print(seq1,seq2,qual1,qual2,offset1,offset2)
##                        mapping[0][9],mapping[0][10] = \
##                                        self.best_seq(seq1,seq2,qual1,qual2,offset1,offset2)

##                    ##sam fields from PE to SE
##                    mapping[0][3] = str(min(int(mapping[0][3]),int(mapping[0][7])))
##                    mapping[0][6],mapping[0][7],mapping[0][8] = "*",str(0),str(0)    
##                    #create new CIGAR string
##                    mapping[0][5] = str(len(mapping[0][9]))+"M"
            
        #Check mismatches
        good_mappings = []
        oligoA_passed_mismatches = 0
        for mapping in mappings:
            seq = mapping[0][9]
            ref = genome[chrom][int(mapping[0][3])-1:\
                                int(mapping[0][3])+len(mapping[0][9])-1].seq
            #creating alignement, and counting matches and mismatches
            alignment,tot_mismatch = self.counting_mismatches(seq,ref)

            #adding tag with mismatching position
            mapping[0] = self.change_add_extra_tag("MD:Z:",alignment,mapping[0])

            #setting allowed mismatched to be depenent from the length of the sequence
            #and given number
            mismatch_limit = max(settings["pseudoSE"]["pseudoSE_allowed_mismatch"],\
#                            ##take minimum from the two:
#                            ###set precentage for the maximum sequenced are
#                            min(settings["pseudoSE"]["pseudoSE_max_read_length"]*2*\
#                                settings["pseudoSE"]["pseudoSE_mismatch_precentage"]/100,\
                        ###precentage from the full length SE sequence
                        len(seq)*settings["pseudoSE"]["pseudoSE_mismatch_precentage"]/100)
            
            #Testing mismatch levels   
            if tot_mismatch > mismatch_limit:
                i = 1
                if int(mapping[0][1]) & 16: #mapped to the negative strand
                    #if oligoA reads are included
                    if settings["pseudoSE"]["pseudoSE_oligoA"]:
                        #testing mismatches without oligoA
                        #and mapping is preserved if mismatches are below threshold
                        while seq[i-1] in {"T","N"}:
                            alignement,tot_mismatch = \
                                            self.counting_mismatches(seq[i+1:],ref[i+1:])
                            i +=1
                            #read considered as not haveing too many mismatches
                            if tot_mismatch <= mismatch_limit:
                                #changing flag just to indicate strand and primary/secondary
                                mapping[0] = self.change_sam_flag(mapping[0],good_mappings)
                                f_oligoA_mismatch.write("\t".join(mapping[0])+"\n")
                                oligoA_mismatch_mappings += 1
                                oligoA_passed_mismatches += 1
                                good_mappings.append(mapping[0])
                                break
                        else: #too many mismatches
                            f_mismatch.write("\t".join(mapping[0])+"\n")
                            mismatched_mappings += 1
                         
                    else: #read considered with too many mismatches
                        f_mismatch.write("\t".join(mapping[0])+"\n")
                        mismatched_mappings += 1

                else: # mapped on positive strand
                     #if oligoA reads are included
                    if settings["pseudoSE"]["pseudoSE_oligoA"]:
                        #testing mismatches without oligoA
                        #and mapping is preserved if mismatches are below threshold
                        while seq[-i] in {"A","N"}:
                            alignement,tot_mismatch = \
                                                self.counting_mismatches(seq[:-i],ref[:-i])
                            i +=1
                            #read considered as not haveing too many mismatches
                            if tot_mismatch <= mismatch_limit:
                                #changing flag just to indicate strand and primary/secondary
                                mapping[0] = self.change_sam_flag(mapping[0],good_mappings)
                                f_oligoA_mismatch.write("\t".join(mapping[0])+"\n")
                                oligoA_mismatch_mappings += 1
                                oligoA_passed_mismatches += 1
                                good_mappings.append(mapping[0])
                                break
                        else: #too many mismatches
                            f_mismatch.write("\t".join(mapping[0])+"\n")
                            mismatched_mappings += 1
                         
                    else: #read considered with too many mismatches
                        f_mismatch.write("\t".join(mapping[0])+"\n")
                        mismatched_mappings += 1

            else: #mismatches below threshold
                #changing flag just to indicate strand and primary/secondary
                mapping[0] = self.change_sam_flag(mapping[0],good_mappings)
                good_mappings.append(mapping[0])

        
        #if all reads are mismatched
        if len(good_mappings) == 0:
            mismatched_reads += 1
            return total_reads,total_mappings,mapped_reads,mappings_counter,\
                many_match_reads,many_match_mappings,mismatched_reads,mismatched_mappings,\
                oligoA_reads,oligoA_mappings,oligoA_mismatch_reads,oligoA_mismatch_mappings,\
                mapping_distribution
        
        #cound oligoA mismatch reads (reads which would fail in oligoA is not considered
        if oligoA_passed_mismatches == len(good_mappings):
            oligoA_mismatch_reads += 1
            
        #Check number of mappings
        ##add tag on match number
        for mapping in good_mappings:
            mapping = self.change_add_extra_tag("NH:i:",len(good_mappings),mapping)

        ##skip read with many mapping and write rest to the file
        if len(good_mappings) > settings["pseudoSE"]["pseudoSE_max_mappings"]:
            for mapping in good_mappings:
                f_many_match.write("\t".join(mapping)+"\n")
            many_match_reads += 1
            many_match_mappings += len(good_mappings)
            return total_reads,total_mappings,mapped_reads,mappings_counter,\
                many_match_reads,many_match_mappings,mismatched_reads,mismatched_mappings,\
                oligoA_reads,oligoA_mappings,oligoA_mismatch_reads,oligoA_mismatch_mappings,\
                mapping_distribution
        else:
            for mapping in good_mappings:
                f_pseudoSE.write("\t".join(mapping)+"\n")
            mapped_reads += 1
            #if len(good_mappings) > 1:
                #print(good_mappings)
            mappings_counter += len(good_mappings)
            print("A",good_mappings)
            print("B",mapping_distribution)
            if len(good_mappings) not in mapping_distribution:
                mapping_distribution[len(good_mappings)] = 1
            else:
                mapping_distribution[len(good_mappings)] += 1
            print("BB",mapping_distribution)

            #saving reads separately which have mismatch A at 3'
            if settings["pseudoSE"]["pseudoSE_oligoA"]:
                oligoA = 0
                for mapping in good_mappings:
                    ref = genome[mapping[2]][int(mapping[3])-1:\
                                             int(mapping[3])-1+\
                                             int(mapping[5].strip("M"))].seq
                    #all reads which do have A-s or N-s in the 3' end and
                    #a single mismatch in this region are considered as
                    #reads with oligoA tail
                    i = 1
                    if int(mapping[1]) & 16:
                        while mapping[9][i-1] in {"T","N"}:
                            if ref[i-1] != "T":
                                f_oligoA.write("\t".join(mapping)+"\n")
                                oligoA_mappings += 1
                                oligoA +=1
                                break
                            i +=1
                    else:
                        while mapping[9][-i] in {"A","N"}:
                            if ref[-i] != "A":
                                f_oligoA.write("\t".join(mapping)+"\n")
                                oligoA_mappings += 1
                                oligoA +=1
                                break
                            i +=1
                #count oligoA reads
                if oligoA == len(good_mappings):
                    oligoA_reads += 1
        return total_reads,total_mappings,mapped_reads,mappings_counter,\
                many_match_reads,many_match_mappings,mismatched_reads,mismatched_mappings,\
                oligoA_reads,oligoA_mappings,oligoA_mismatch_reads,oligoA_mismatch_mappings,\
                mapping_distribution

    def convert_to_phred33(self,phred_string):
        '''
        Converts phred+64 score to phred+64 score
        '''
        phred33 = "$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"
        phred64 = 'CDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghi'
        phred_dict = { phred64[i]:phred33[i] for i in range(78) if i < 39}
        return "".join([phred_dict[base] for base in reversed(phred_string)])
        
    def best_seq(self,seq1,seq2,qual1,qual2,offset1,offset2):
        '''
        Combines two overlaping sequencies according to score.
        Gives also combined score.
        '''
        new_seq = ""
        new_qual = ""

        if offset1 >= 0: 
            for x in range (0,abs(offset1)): #5' no overlap region
                new_seq += seq2[x]
                new_qual += qual2[x]
            for x in range (0,len(seq1)-max(offset2,0)): #overlap region
                best =  self.best_base(seq1,seq2,qual1,qual2,x,offset1,offset2)
                new_seq += best[0]
                new_qual += best[1]
            if offset2 > 0: #3' no overlap region
                new_seq += seq1[-offset2:]
                new_qual += qual1[-offset2:]            
            elif offset2 < 0:
                new_seq += seq2[offset2:]
                new_qual += qual2[offset2:]
        elif offset1 < 0:
            for x in range (0,abs(offset1)): #5' no overlap region
                new_seq += seq1[x]
                new_qual += qual1[x]
            for x in range (0,len(seq2)+min(offset2,0)): #overlap region
                best =  self.best_base(seq2,seq1,qual2,qual1,x,-offset1,offset2)
                new_seq += best[0]
                new_qual += best[1]
            if offset2 > 0: #3' no overlap region
                new_seq += seq1[-offset2:]
                new_qual += qual1[-offset2:]            
            elif offset2 < 0:
                new_seq += seq2[offset2:]
                new_qual += qual2[offset2:]
        return new_seq, new_qual

    def best_base(self,seq1,seq2,qual1,qual2,pos,offset1,offset2):
        '''
        Identifies nucleotide with better score.
        '''
        phred33="!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"
        pos_qual1=qual1[pos]
        pos_qual2=qual2[pos+offset1]
        phred33_pos_qual1 = phred33.find(pos_qual1)
        phred33_pos_qual2 = phred33.find(pos_qual2)
        if phred33_pos_qual1 >= phred33_pos_qual2:
            return(seq1[pos], pos_qual1)
        else:
            return(seq2[pos+offset1], pos_qual2)
                                                               
    def change_sam_flag(self,mapping,good_mappings):
        '''
        Changing FLAG to indicate only strand and primary or secondary alignment
        '''
        if len(good_mappings) == 0: #primaryalignment
            if int(mapping[1]) & 16:
                mapping[1] = str(16)
            else:
                mapping[1] = str(0)
                
        else: #secondary alignment
            if int(mapping[1]) & 16:
                mapping[1] = str(16+256)
            else:
                mapping[1] = str(0+256)
        return mapping
    
    def counting_mismatches(self,seq,ref):
        '''
        Counts mismatches between to input sequencies
        and gives alingment as an MD tag 
        '''
        alignement = ""
        match = 0
        tot_mismatch = 0
        for x in range(0,len(ref)):
            if ref[x] == seq[x]:
                match += 1
                if x == len(seq)-1:
                    alignement += str(match)
            #now considering "N" as match but also contributing to mismatches
            elif seq[x] == "N": 
                match +=1
                tot_mismatch += 0.25
                if x == len(seq)-1:
                    alignement += str(match)
            else:
                tot_mismatch +=1
                alignement += str(match)+ref[x]
                match = 0 
        return alignement,tot_mismatch

    def change_add_extra_tag(self,tag,tag_value,list_sam):
        '''
        Change or add (if missing) tag to read in sam format    
        '''
        extra_tags = list_sam[11:]
        tag_index = [i for i, j in enumerate(extra_tags) if j.startswith(tag)]
        #tag does not exist
        if tag_index == []:
            list_sam.append(tag+str(tag_value))
        else:
            for i,x in enumerate(tag_index):
                extra_tags[tag_index[i]] = tag+str(tag_value)
                list_sam = list_sam[:11]+extra_tags
        return list_sam

    def write_statistics(self,settings,library,total_reads,mismatched_reads,many_match_reads,
                            mapped_reads,oligoA_reads,oligoA_mismatch_reads,total_mappings,
                            mismatched_mappings,many_match_mappings,mappings_counter,
                            oligoA_mappings,oligoA_mismatch_mappings,mapping_distribution):
        '''
        Writes data about pseudoSE generation.
        '''
        info_file = os.path.join(settings["--output"],"pseudoSE","pseudoSE_info",\
                                  library+"_pseudoSEinfo.log")
        f_info = open(info_file,"w")
        
        f_info.write("Input reads:\t"+str(total_reads) + "\n")
        f_info.write("Mismatched reads:\t"+str(mismatched_reads) + "\n")
        f_info.write("Reads with too many matches:\t"+str(many_match_reads) + "\n")
        f_info.write("Passed reads:\t"+str(mapped_reads) + "\n")
        if settings["pseudoSE"]["pseudoSE_oligoA"]:
            f_info.write("Passed oligoA reads:\t"+str(oligoA_reads) + "\n")
            f_info.write("Passed oligoA reads with mismatches:\t"+str(oligoA_mismatch_reads) + \
                         "\n\n")
        else:
            f_info.write("\n")
        f_info.write("Input genomic matches:\t"+str(total_mappings) + "\n")
        f_info.write("Mismatched genomic matches:\t"+str(mismatched_mappings) + "\n")
        f_info.write("Genomic matches of reads with too many matches:\t"+\
                     str(many_match_mappings) + "\n")
        f_info.write("Passed genomic matches:\t"+str(mappings_counter) + "\n")
        if settings["pseudoSE"]["pseudoSE_oligoA"]:
            f_info.write("Passed oligoA genomic matches:\t"+str(oligoA_mappings) + "\n")
            f_info.write("Passed oligoA genomic matches with mismatches:\t"+\
                         str(oligoA_mismatch_mappings) + "\n\n")
        else:
            f_info.write("\n")

        f_info.write("Mapping distribution:" + "\n")
        for element in sorted(mapping_distribution):
            f_info.write(str(element) + "\t" + str(mapping_distribution[element]) + "\n")
        if many_match_reads > 0:
            f_info.write(">"+str(settings["pseudoSE"]["pseudoSE_max_mappings"]) + "\t" +\
                         str(many_match_reads) + "\n")
        f_info.close()

    def test_errors(self,settings):
        '''
        Testing existence and errors in info files
        '''
        info_folder = os.path.join(settings["--output"],"pseudoSE","pseudoSE_info")
        if settings["pseudoSE"]["pseudoSE_oligoA"]:
            wanted_read_beginnings = {"Input","Passed","Mismatched","Reads","Genomic"}
        else:
            wanted_read_beginnings = {"Input","Passed","Mismatched","Reads","Genomic"}
        if len(os.listdir(info_folder)) == 0:
            sys.exit('Task "pseudoSE" incomplete! Infofiles missing in folder ' + info_folder)
        for file in os.listdir(info_folder):
            with open(os.path.join(info_folder,file)) as f_in:
                read_beginnings = set()
                for line in f_in:
                    read_beginnings.add(line.split(" ")[0].split("\t")[0].strip())
                if not wanted_read_beginnings.issubset(read_beginnings):
                    sys.exit('Error in task "pseudoSE", check infofile: '+\
                             os.path.join(info_folder,file))
