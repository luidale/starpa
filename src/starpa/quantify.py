#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Count reads for clustered pp-s.
Create annotation table if missing.
Count reads by annotation.
Count reads for different annotation groups.
Normalize counts.
Collect data.

To do:
1)Consensus sequence: double degenerated nucleotides
2)Create contig html files
3)Is the splitting of the strand nessecary?

'''

import os
import bisect
import copy
import sys
import multiprocessing as mp
from collections import defaultdict

from starpa.annotation import annotation 

from pyfaidx import Fasta

class quantify():
    def __init__(self,settings,first_task):
        self.settings = settings
        self.make_folder(settings)
        self.quantify(settings,first_task)
        self.test_errors(settings)

    def make_folder(self,settings):
        '''
        Make output folder for task if needed
        '''
        if not os.path.exists(os.path.join(settings["--output"],"quantify")):
            os.makedirs(os.path.join(settings["--output"],"quantify"))
        if not os.path.exists(os.path.join(settings["--output"],"quantify","libraries")):
            os.makedirs(os.path.join(settings["--output"],"quantify","libraries"))
        for library in sorted(settings["libraries"]):
            if not os.path.exists(os.path.join(settings["--output"],"quantify","libraries",library)):
                os.makedirs(os.path.join(settings["--output"],"quantify","libraries",library))
        if not os.path.exists(os.path.join(settings["--output"],"quantify","selected_pps")):
            os.makedirs(os.path.join(settings["--output"],"quantify","selected_pps"))
        if not os.path.exists(os.path.join(settings["--output"],"quantify","collected_statistics")):
            os.makedirs(os.path.join(settings["--output"],"quantify","collected_statistics"))
            
    def quantify(self,settings,first_task):
        '''
        Quantifies pp-s and collects data
        '''
        overlap = settings["overlap_range"]
        size_range = settings["size_range"]
        
        #CHECK ANNOTATION FILE
        #is it BED file
        if settings["quantify"]["quantify_annotation_file"][-3:] != "BED":
            #is converted BED file existing
            if not os.path.isfile(settings["quantify"]["quantify_annotation_file"][:-3]+"BED"):
                if not os.path.isfile(settings["quantify"]["quantify_annotation_file"][:-4]+"BED"):
                    #convert GFF or GFF3 to BED format
                    annotation(settings)
        
        annotation_BED_file = ".".join(settings["quantify"]["quantify_annotation_file"].\
                                       split(".")[:-1])+".BED"            
        settings["quantify"]["quantify_annotation_file"] = annotation_BED_file

        #GET PP BIOTYPES
        self.get_pp_biotype(settings,first_task)

        #GET LISTS FOR ANNOTATION ELEMENTS AND GROUPS
        #create lst of genes with their positions
        print("\tCreate a list for annotation elements")
        gene_list = self.genes_from_BED(settings["quantify"]["quantify_annotation_file"])      
        #Version 4 -creates lists with start and end positions of genes, and names of the genes
        #at those positions
        gene_pos_list, gene_name_list = self.genes_to_positions4(gene_list)
        
        #SPLIT FILE WITH PP-s
        strand_list = ["For","Rev"]
        print("\tSplit pp-file")
        self.split_by_strand_BED(settings)
##        for strand_name in strand_list:
##            pp_stranded_file = os.path.join(settings["--output"],"quantify",\
##                                              "pp_clustered"+"_"+strand_name+".BED")
##            self.fragment_BED(pp_stranded_file,overlap,size_range)

        #GET LISTS FOR PP-s
        #create lst of genes with their positions
        print("\tCreate list for pps")
        pp_list = {}
        for strand_name in strand_list:
            if strand_name == "For":
                strand = "+"
            else:
                strand = "-"
            pp_list[strand] = self.genes_from_BED(os.path.join(settings["--output"],\
                                                                    "quantify",\
                                              "pp_clustered"+"_"+strand_name+".BED"))
                    
        #Version 4 -creates lists with start and end positions of genes, and names of the genes
        #at those positions
        pp_pos_list = {}
        pp_name_list = {}
        for strand_name in strand_list:
            if strand_name == "For":
                strand = "+"
            else:
                strand = "-"
            pp_pos_list[strand], pp_name_list[strand] = \
                                      self.genes_to_positions4(pp_list[strand])
            
        #GET DICTIONARY FOR PP-S
        print("\tMake pp dictionary")
        pp_data = self.make_pp_dictionary(settings,first_task)

        #ADD GENOMIC SEQUENCE TO PP DICT
        pp_data = self.add_genomic_sequence(settings,pp_data)

        #COUNT READS FOR PP-s
        ##libraries are proccessed in paralel
        print("\tCount reads for pp-s")
        pool = mp.Pool(processes=settings["CPUs"])
        results = {library:pool.apply_async(self.parse_reads, \
                          args = (settings,library,strand_list,first_task,pp_data,
                                  gene_list,gene_pos_list,gene_name_list,
                                  pp_list,pp_pos_list,pp_name_list)) \
                          for library in sorted(settings["libraries"])}
        pp_dict = {library:p.get() for library,p in results.items()}
        pool.close()
        pool.join()

        #REMOVE SPLITTED PP_FILE
        for strand_name in strand_list:
            os.remove(os.path.join(settings["--output"],"quantify",\
                                   "pp_clustered"+"_"+strand_name+".BED"))

        #WRITE COMBINED BIOTYPE COUNT DATA
        ##get biotypes
        biotypes = [x for x in sorted(pp_dict[list(settings["libraries"])[0]]["biotype_counts"]) \
                                                                            if x != ("total")]
        ##write data
        with open(os.path.join(settings["--output"],"quantify",\
                               "collected.biotype_annotation.statistics"),"w") as f_out:
            f_out.write("\t".join(["#libraries"]+["total"]+biotypes)+"\n")
            for library in sorted(settings["libraries"]):
                f_out.write(library)
                f_out.write("\t"+str(round(pp_dict[library]["biotype_counts"]["total"],1)))
                for biotype in biotypes:
                    f_out.write("\t"+str(round(sum(float(x) for x in \
                                            pp_dict[library]["biotype_counts"][biotype]),1)))
                f_out.write("\n")
                
        #Write counts
        print("\tWrite normalized pp counts")
        comment = "Read counts given in absolute numbers"
        self.write_pp_recalc_file(settings,pp_dict,"total",comment)
        comment = "Read counts given in RPM-s"
        self.write_pp_recalc_file(settings,pp_dict,"RPM",comment)
        comment = "Read counts given in RPM-s of biotype, all antisense pp-s are considered as "+\
                                                                                      "intergenic"
        self.write_pp_recalc_file(settings,pp_dict,"biotype_RPM",comment)
        comment = "Read counts given in RPM-s of groupped biotypes, non groupped are "+\
                          ",".join(settings["quantify"]["quantify_non_groupped_biotypes"])
        self.write_pp_recalc_file(settings,pp_dict,"groupped_biotype_RPM",comment)

        #Write data
        print("\tWrite pp data")
        comment = "Consensus sequence +/- "+str(settings["non_overlap"])+" nucleotides "+\
                                "(extra nucleotides separated by *): eg. CT*CTAGTCAGTC*T "                
        self.write_pp_data_file(settings,pp_dict,"consens_seq",comment)
        comment = "Comparative genominc sequence ('.' when match with consensus sequence) "+\
                        "+/- "+str(settings["non_overlap"])+\
                        " nucleotides (extra nucleotides separated by *): eg. ..*...C.....*A "                
        self.write_pp_data_file(settings,pp_dict,"genomic_seq",comment)
        comment = "Quality of consensus sequence +/- "+str(settings["non_overlap"])+\
                  " nucleotides (extra nucleotides separated by *): eg. 89*999999999*8 " 
        self.write_pp_data_file(settings,pp_dict,"consens_qual",comment)
        comment = "Coverage of reads of pp +/- "+str(settings["non_overlap"])+" nucleotides "+\
                                 "(extra nucleotides separated by *): eg. 0;2;*;50;50...50...50;49;*10;0"
        self.write_pp_data_file(settings,pp_dict,"coverage",comment)
        comment = "Genomic uniqness of pp calculated by taking an average of mapping number of "+\
                                                                                "all reads of the pp"
        self.write_pp_data_file(settings,pp_dict,"uniqness",comment)
        comment = "Relative coverage of calculated by calculating coverage of pp divided by "+\
                                                                             "coverage in the region"
        self.write_pp_data_file(settings,pp_dict,"rel_cov",comment)

        #Write pp-s selected by counts
        comment = "Read counts given in absolute numbers"
        self.write_pp_recalc_file_selected(settings,pp_dict,"total",comment)
        comment = "Read counts given in RPM-s"
        self.write_pp_recalc_file_selected(settings,pp_dict,"RPM",comment)

        #Collect statistics
        self.collect_statistics(settings)
        
    def collect_statistics(self,settings):
        '''
        Collects library specific data to statistics folder
        '''
        tasks = ["trim", "align", "sam_sort", "pseudoSE", "identify"]
        folder = []
        folder.append(os.path.normpath(os.path.join(setting["--input"],"..")))
        folder.append(os.path.normpath(settings["--output"]))
        for libraries in settings["libraries"]:
            with open(os.path.join(settings["--output"],"quantify","collected_statistics",\
                                   "collected_stat_"+library+".log"),"w") as f_out:
                for task in tasks:
                    for folder_path in folder:
                        if os.path.isdir(os.path.join(folder_path,path)):
                            f_out.write("#"+task+"\n")
                            if task == "identify":
                                for strand in ["For","Rev"]:
                                    f_out.write("#"+strand+"\n")
                                    with open(os.path.join(folder_path,path,task+"_info",
                                                           library+"_"+strand+"_"+task+"info.log")) as f_in:
                                        for line in f_in:
                                            f_out.write(line)
                                        f_out.write("\n")                                
                            else:
                                with open(os.path.join(folder_path,path,task+"_info",
                                                       library+"_"+task+"info.log")) as f_in:
                                    for line in f_in:
                                        f_out.write(line)
                                    f_out.write("\n")
                            break

    def genes_from_BED(self,annotation_file):
        '''
        Creats dictionary for chromosomes containing lists with gene information. Info comes
        from BED so it is 0 based end exclusiv.
        '''
        #print("\tCreate gene list")
        gene_list = {}
        f1 = open(annotation_file)
        i = 1
        for line in f1:
            #coordination system 0-based, end exclusive
            line = line.strip().split("\t")
            chrom,start,end,gene,strand,biotype = line[0],int(line[1]),int(line[2]),line[3],line[5],line[6]
            if chrom not in gene_list:
                gene_list[chrom] = [[0,0,"","","intergenic",0,0,0,0]]
            #[...sense exon,antisense exon, undirected region, gene_id]:
            gene_list[chrom].append([start,end,strand,gene,biotype,0,0,0,i]) 
            i += 1
            #print(gene_list[chrom][-1])
        f1.close()  
        return gene_list

    def genes_to_positions4(self,gene_list):
        '''
        Creates dictionary with positions of the genes refering to the name of the genes.
        Input: 0-based end exclusive
        Output: 0-based and end inclusive
        '''
        #print("\tAdding genes to positions")
        gene_dict = {}
        single_pos_element = set()
        for chrom in gene_list:
            '''creates dictionary with start and end positions of the elements
            {1:{1}, 10:{1}, 15:{2}, 20:{2}, 8:{3}, 12:{3}}
            '''    
            for gene in gene_list[chrom][1:]:
                if chrom not in gene_dict:
                    gene_dict[chrom] = {}
                if int(gene[0]) not in gene_dict[chrom]:
                    gene_dict[chrom][int(gene[0])] = {gene[8]}
                else:
                    gene_dict[chrom][int(gene[0])].add(gene[8])
                #single position elements
                if int(gene[1])-int(gene[0]) == 1:
                    single_pos_element.add(gene[8])
                    
                if int(gene[1])-1 not in gene_dict[chrom]:
                    gene_dict[chrom][int(gene[1])-1] = {gene[8]}
                else:
                    gene_dict[chrom][int(gene[1])-1].add(gene[8]) 

            '''correct overlaping elements
            #if elements overlap and the one elements starts in the middle of the other
            #then the elements get cross labeled in the positions
            {1:{1}, 8:{1,3}, 10:{1,3}, 12:{3}, 15:{2}, 20:{2}}
            '''
            genes_started = set()
            for x in sorted(gene_dict[chrom]):
                #if no started elements
                if genes_started == set():
                    genes_started = copy.deepcopy(gene_dict[chrom][x])
                #if some started elements
                else:
                    #nothing ending - nothing common between started and ongoing
                    if genes_started.intersection(gene_dict[chrom][x]) == set():
                        #adding started
                        genes_started = genes_started.union(gene_dict[chrom][x])
                        gene_dict[chrom][x] = copy.deepcopy(genes_started)
                    #something ending
                    else:
                        #for each ending
                        old_genes_started = copy.deepcopy(genes_started)
                        genes_started = genes_started.union(gene_dict[chrom][x])
                        #removing from started list those which are common for position and started
                        for z in old_genes_started.intersection(gene_dict[chrom][x]):
                            #print(z)
                            genes_started.remove(z)    
                        gene_dict[chrom][x] = genes_started.union(gene_dict[chrom][x])
                #single pos elements are removed from starting elements,
                #otherwise their end are searched
                if single_pos_element.intersection(genes_started) != set():
                    for y in single_pos_element.intersection(genes_started):
                        genes_started.remove(y)

        #convert dictionary to lists
        gene_pos_list = {}
        gene_name_list = {}
        for chrom in gene_dict:
            gene_pos_list[chrom] = [0]
            gene_name_list[chrom] = [set()]
            for pos in sorted(gene_dict[chrom]):
                gene_pos_list[chrom].append(pos)
                gene_name_list[chrom].append(gene_dict[chrom][pos])
            gene_name_list[chrom].append(set())
        return gene_pos_list, gene_name_list

    def split_by_strand_BED(self,settings):
        '''
        Split BED by strand
        '''
        input_file = os.path.join(settings["--output"],"quantify","pp_clustered_biotype.BED")
        output_file = os.path.join(settings["--output"],"quantify","pp_clustered")
        f_in = open(input_file)
        f_out_pos = open(output_file+"_For.BED","w")
        f_out_neg = open(output_file+"_Rev.BED","w")
        for line in f_in:
            while line[0] == "#": #skip header
                line = f_in.realine()
            if line.strip().split("\t")[5] == "+":
                f_out_pos.write(line)
            elif line.strip().split("\t")[5] == "-":
                f_out_neg.write(line)
        f_out_pos.close()
        f_out_neg.close()
        f_in.close()

    def get_pp_biotype(self,settings,first_task):
        '''
        Get biotype of pp-s
        '''
        #match pp-s and genes
        print("\tMatch pp-s and genes")
        if first_task == "quantify":
            pp_input_file = os.path.join(settings["--input"],"pp_clustered.BED")
        else:
            pp_input_file = os.path.join(settings["--output"],"cluster","pp_clustered.BED")
            
        pp_gene_match_command = (
                        settings["bedtools_call"], "intersect",
                        "-a", pp_input_file,
                        "-b", ".".join(settings["quantify"]["quantify_annotation_file"]\
                                       .split(".")[:-1])+".BED",
                        "-wao",
                        ">", os.path.join(settings["--output"],"quantify",\
                                          "pp_clustered_biotype_match.BED")
                        )
        os.system(" ".join(pp_gene_match_command))

        #read in pp_gene and gene matches
        pp_match_biotype_dic = {}
        pp_match_gene_dic = {}
        with open(os.path.join(settings["--output"],"quantify",\
                                          "pp_clustered_biotype_match.BED")) as f_pp_match:
            for line in f_pp_match:
                line = line.strip().split("\t")
                pp_strand,pp_name,gene_strand = line[5],line[3],line[11]
                biotype,gene = line[12],":".join([line[9],line[12],line[13]])
                #test strand
                if pp_strand not in pp_match_biotype_dic:
                    pp_match_biotype_dic[pp_strand] = {}
                if pp_strand not in pp_match_gene_dic:
                    pp_match_gene_dic[pp_strand] = {}
                #test pp
                if pp_name not in pp_match_biotype_dic[pp_strand]:
                    pp_match_biotype_dic[pp_strand][pp_name] = {}
                    pp_match_biotype_dic[pp_strand][pp_name]["sense"] = []
                    pp_match_biotype_dic[pp_strand][pp_name]["antisense"] = []
                if pp_name not in pp_match_gene_dic[pp_strand]:
                    pp_match_gene_dic[pp_strand][pp_name] = {}
                    pp_match_gene_dic[pp_strand][pp_name]["sense"] = []
                    pp_match_gene_dic[pp_strand][pp_name]["antisense"] = []
                #test was there a match
                if line[6] == ".":
                    continue
                #select proper orientation of the gene
                if gene_strand == pp_strand:
                    pp_match_biotype_dic[pp_strand][pp_name]["sense"].append(biotype)
                    pp_match_gene_dic[pp_strand][pp_name]["sense"].append(gene)
                else:
                    pp_match_biotype_dic[pp_strand][pp_name]["antisense"].append(biotype)
                    pp_match_gene_dic[pp_strand][pp_name]["antisense"].append(gene)
        #os.remove(os.path.join(input_folder, pp_file)[:-4]+"_biotype.BED")

        ##give each pp a biotype and gene
        pp_dic = {}
        for pp_strand in pp_match_biotype_dic:
            for pp_name in pp_match_biotype_dic[pp_strand]:
                pp_dic[pp_name] = {}
                pp_dic[pp_name]["biotype"] = set()
                pp_dic[pp_name]["gene"] = set()

                #testing biotypes in sense strand
                if pp_match_biotype_dic[pp_strand][pp_name]["sense"] != []:
                    for biotype in pp_match_biotype_dic[pp_strand][pp_name]["sense"]:
                        pp_dic[pp_name]["biotype"].add(biotype)
                    for gene in pp_match_gene_dic[pp_strand][pp_name]["sense"]:    
                        pp_dic[pp_name]["gene"].add(gene)
                else: #if nothing in sense strand, test antisense strand
                    for biotype in pp_match_biotype_dic[pp_strand][pp_name]["antisense"]:
                        pp_dic[pp_name]["biotype"].add("antisense_"+biotype)
                    for gene in pp_match_gene_dic[pp_strand][pp_name]["antisense"]:    
                        pp_dic[pp_name]["gene"].add("antisense:"+gene)                        
                #if no biotype given name as intergenic
                if pp_dic[pp_name]["biotype"] == set():
                    pp_dic[pp_name]["biotype"].add("intergenic")
                    pp_dic[pp_name]["gene"].add("intergenic") 
                pp_dic[pp_name]["biotype"] = sorted(list(pp_dic[pp_name]["biotype"]))
                pp_dic[pp_name]["gene"] = sorted(list(pp_dic[pp_name]["gene"]))
                
        ##write biotype in to the file
        print("\tWrite pp biotype file") 
        pp_bio_out = open(os.path.join(settings["--output"],"quantify",\
                                          "pp_clustered_biotype.BED"),"w")
        with open(pp_input_file) as pp_in:
            for line in pp_in:
                line = line.strip().split("\t")
                line.append(";".join(pp_dic[line[3]]["biotype"]))
                pp_bio_out.write("\t".join(line)+"\n")
        pp_bio_out.close()
        print("\tWrite pp gene file") 
        pp_bio_out = open(os.path.join(settings["--output"],"quantify",\
                                          "pp_clustered_gene.BED"),"w")
        with open(pp_input_file) as pp_in:
            for line in pp_in:
                line = line.strip().split("\t")
                line.append(";".join(pp_dic[line[3]]["gene"]))
                pp_bio_out.write("\t".join(line)+"\n")
        pp_bio_out.close()


    def make_pp_dictionary(self,settings,first_task):
        '''
        Makes empty dictionary to save data for pps.
        '''
        #file name
        input_file_biotype = os.path.join(settings["--output"],"quantify","pp_clustered_biotype.BED")
        input_file_gene = os.path.join(settings["--output"],"quantify","pp_clustered_gene.BED")        
        with open(input_file_biotype) as f_in:
            pp_data = {}
            for line in f_in:
                line = line.strip().split("\t")
                pp_data[line[3]] = {}
                #chrom,start,end,strand,biotype:
                pp_data[line[3]]["data"] = [line[0],int(line[1]),int(line[2]),line[5]]
                pp_data[line[3]]["biotype"] = line[6].split(";")
                pp_data[line[3]]["count"] = {}
                pp_data[line[3]]["count"]["total"] = [0]
                pp_data[line[3]]["uniqness"] = 0
                #[A,C,G,T,N] for each position:
                pp_data[line[3]]["sequencies"] = [[0 for x in range(5)] for x in \
                                                  range(int(line[2])-int(line[1])+2*settings["non_overlap"])]
                pp_data[line[3]]["coverage"] = [0 for x in range(int(line[2])-int(line[1])+2*settings["non_overlap"])]
                pp_data[line[3]]["consens_seq"] = ""
                pp_data[line[3]]["consens_qual"] = ""
                pp_data[line[3]]["rel_cov"] = ""
                pp_data[line[3]]["total_cov"] = 0
                pp_data[line[3]]["genomic_seq"] = ""
                pp_data[line[3]]["raw_genomic_sequence"] = ""
        #add gene data
        with open(input_file_gene) as f_in:
            for line in f_in:
                line = line.strip().split("\t")
                pp_data[line[3]]["gene"] = line[6].split(";")
        return pp_data

    def add_genomic_sequence(self,settings,pp_data):
        '''
        Adds genomic sequence for each pps
        '''
        genome = Fasta(settings["genome"],one_based_attributes=False)
        for pp_name in pp_data:
            chrom,start,end,strand = (x for x in pp_data[pp_name]["data"][:4])
            pp_gen_seq = genome[chrom][int(start)-settings["non_overlap"]:\
                                        int(end)+settings["non_overlap"]].seq
            pp_data[pp_name]["raw_genomic_sequence"] = pp_gen_seq
        return pp_data
    
    def parse_reads(self,settings,library,strand_list,first_task,pp_data,
                                  gene_list,gene_pos_list,gene_name_list,
                                  pp_list,pp_pos_list,pp_name_list):
        '''
        Parsing reads and counting them for pp-s and genes
        '''
        print("\t\t"+library,"Parsing mappings")        
        gene_list_lib = copy.deepcopy(gene_list)
        #set file names
        if first_task in {"cluster","quantify"}:
            input_sam = os.path.normpath(os.path.join(settings["--input"],\
                                     settings["quantify"]["quantify_sam_file_location"],\
                                 library+settings["quantify"]["quantify_sam_file_suffix"]))
        elif first_task == "identify":
            input_sam = os.path.normpath(os.path.join(settings["--input"],\
                                 library+settings["quantify"]["quantify_sam_file_suffix"]))            
        else:
            input_sam = os.path.join(settings["--output"],"pseudoSE",\
                                 library+"_pseudoSE.sam")

        #parse file   
        with open(input_sam) as f_in:
            for line in f_in:
                #skip header
                while line.startswith("@"):
                    line = f_in.readline()
                #set strand
                line_split = line.strip().split("\t")
                if int(line_split[1]) & 16:
                    strand = "-"
                else:
                    strand = "+"
                #SAM coordination system (1-based) converted to 0-based end inclusive
                start = int(line_split[3])-1
                end = start + abs(int(line_split[5].strip("M")))-1
                chrom = line_split[2]
                #matching reads to genes
                gene_list_lib = self.reads_to_genes4(line,gene_list_lib,gene_pos_list,\
                                                     gene_name_list,strand,start,end,chrom)
                #matching reads to pps
                pp_data = self.reads_to_pps4(settings,line_split,pp_list[strand],pp_data,\
                                                 pp_pos_list[strand],pp_name_list[strand],\
                                                 start,end,chrom)

        #calculate coverage, uniqness, consensus sequence, consensus quality,genomic sequence
        print("\t\t"+library,"Calculate pp properties")    
        for pp_name in pp_data:
            pp_data[pp_name] = self.calculate_pp_data(settings,pp_data[pp_name])

        #write pp data to the file
        print("\t\t"+library,"Write data")
        with open(os.path.join(settings["--output"],"quantify","libraries",library,\
                               "pp_clustered_"+library+"_counted.BED"),"w") as f_out:
            f_out.write("\t".join(["#Chromosome","Start","End","PP_name","Score","Strand",\
                                   "Count","Consensus_seq","Consensus_qual","Genomic_qual",\
                                   "Coverage","Uniqness","Relative_coverage"])+"\n")
            for pp_name in sorted(pp_data,key=lambda y:(pp_data[y]["data"])):
                #pp-s with count
                if pp_data[pp_name]["count"]["total"][0] != 0:
                    #[chrom,start,end,name,score,strand,
                    # count,consensus_sequence,consensus_quality,coverage,
                    # uniqness,relative_coverage]# coverage as 1;15;17;20 etc
                    f_out.write("\t".join([str(x) for x in pp_data[pp_name]["data"][:3]]+\
                                          [pp_name,"1",pp_data[pp_name]["data"][3],\
                                           str(pp_data[pp_name]["count"]["total"][0]),\
                                           pp_data[pp_name]["consens_seq"],\
                                           pp_data[pp_name]["consens_qual"],\
                                           pp_data[pp_name]["genomic_seq"],\
                                           pp_data[pp_name]["coverage"],\
                                           str(pp_data[pp_name]["uniqness"]),\
                                           str(pp_data[pp_name]["rel_cov"])])+"\n")
                else:
                    f_out.write("\t".join([str(x) for x in pp_data[pp_name]["data"][:3]]+\
                                          [pp_name,"1",pp_data[pp_name]["data"][3],\
                                           str(pp_data[pp_name]["count"]["total"][0]),\
                                           ".",".",".",".",".","."])+"\n")
        #Write annotation statistics
        print("\t\t"+library,"Write data")
        self.write_annotation_statistics(gene_list_lib,\
                                    os.path.join(settings["--output"],"quantify","libraries",\
                                                 library,library+".gene_annotation.statistics"),
                                    os.path.join(settings["--output"],"quantify","libraries",\
                                                 library,library+".biotype_annotation.statistics"))
        #Get biotype counts
        biotype_counts = {}
        biotype_strand_counts = {}
        for chrom in gene_list_lib:
            for gene in gene_list_lib[chrom]:
                #leaves out the number which was gene order number:
                if gene[4] not in biotype_counts:
                    biotype_counts[gene[4]] = gene[5]
                    biotype_strand_counts[gene[4]] = gene[5:-1]
                else:
                    biotype_counts[gene[4]] += gene[5]
                    biotype_strand_counts[gene[4]] = [biotype_strand_counts[gene[4]][i] + x for i,x \
                                                      in enumerate(gene[5:-1])]
                if "intergenic" not in biotype_counts:
                    biotype_counts["intergenic"] = sum(gene[6:-1])
                else:
                    biotype_counts["intergenic"] += sum(gene[6:-1])
        biotype_counts["total"] = sum([biotype_counts[x] for x in biotype_counts])
        biotype_strand_counts["total"] = sum([sum(biotype_strand_counts[x]) for x \
                                              in biotype_strand_counts])
        pp_data["biotype_counts"] = biotype_strand_counts
        
        #Recalculate pp count data
        #recalculate pp count values
        #convert read numbers to RPM
        print("\t\t"+library,"Convert pp-counts")#
        for pp_name in pp_data:
            if pp_name == "biotype_counts":
                continue
            pp_data[pp_name]["count"]["RPM"] = {}
            pp_data[pp_name]["count"]["RPM"] = \
                        [str(round(pp_data[pp_name]["count"]["total"][0]/\
                                   biotype_counts["total"]*1000000,5))]

        #convert read numbers to RPM of biotype
        for pp_name in pp_data:
            if pp_name == "biotype_counts":
                continue
            biotypes = pp_data[pp_name]["biotype"]
            pp_data[pp_name]["count"]["biotype_RPM"] = []
            for biotype in biotypes:
                #sense and antisense repeat_regions will be repeat_regions
                if biotype.find("repeat_region") != -1: 
                    biotype = "repeat_region"
                #considers all antisense biotypes as intergenic
                if biotype.find("antisense") != -1: 
                    biotype = "intergenic"
                if float(biotype_counts[biotype]) == 0:
                    pp_data[pp_name]["count"]["biotype_RPM"].append(str(0))
                else:
                    try:
                        pp_data[pp_name]["count"]["biotype_RPM"].append(\
                            str(round(int(pp_data[pp_name]["count"]["total"][0])/\
                                      float(biotype_counts[biotype])*1000000,5)))
                    except:#for testing
                        print(library, biotype,biotype_counts[biotype])

        #convert read numbers to RPM of groupped biotype:
        #1) tRNA
        #2) rRNA
        #3) repeat regions
        #4) all others
        ##get read count for groupped biotypes 
        groupped_read_count = self.get_grupped_biotype_reads(\
                            biotype_counts,settings["quantify"]["quantify_non_groupped_biotypes"])
        for pp_name in pp_data:
            if pp_name == "biotype_counts":
                continue
            pp_data[pp_name]["count"]["groupped_biotype_RPM"]=[]
            biotypes = pp_data[pp_name]["biotype"]
            for biotype in biotypes:
                #sense and antisense repeat_regions will be repeat_regions
                if biotype.find("repeat_region") != -1: 
                    biotype = "repeat_region"
                #considers all antisense biotypes as intergenic
                if biotype.find("antisense") != -1: 
                    biotype = "intergenic"
                if biotype in settings["quantify"]["quantify_non_groupped_biotypes"]:
                    if float(biotype_counts[biotype]) == 0:
                        pp_data[pp_name]["count"]["groupped_biotype_RPM"].append(str(0))
                    else:
                        try:
                            pp_data[pp_name]["count"]["groupped_biotype_RPM"].append(\
                                str(round(int(pp_data[pp_name]["count"]["total"][0])/\
                                          float(biotype_counts[biotype])*1000000,5)))
                        except:
                            print(library, biotype,biotype_counts[biotype])
                else: #all other biotypes are handled as groupped
                    if float(biotype_counts[biotype]) == 0:
                        pp_data[pp_name]["count"]["groupped_biotype_RPM"].append(str(0))
                    else:
                        try:
                            pp_data[pp_name]["count"]["groupped_biotype_RPM"].append(\
                                str(round(int(pp_data[pp_name]["count"]["total"][0])/\
                                          groupped_read_count*1000000,5)))
                        except:
                            print(library, biotype,groupped_read_count)                
        
        print("\t\t"+library,"Done")  
        return pp_data

    def reads_to_genes4(self,mapping,gene_list,gene_pos_list,gene_name_list,strand,start,end,chrom):
        '''
        Looks up genes overlaping with mapping and puts mapping info to gene_list
        '''
        
        #get the annotation elements in the positions of the mapping
        ##get positions which are also positions of the annotation elements
        ##and gets directly gene names        
        index1=bisect.bisect_left(gene_pos_list[chrom],start)
        index2=bisect.bisect(gene_pos_list[chrom],end)
        #print("b", index1,index2)
        if index1 == index2:
            genes_on_mapping = gene_name_list[chrom][index1].\
                               intersection(gene_name_list[chrom][index1-1])
        else:
            genes_on_mapping = set.union(*map(set,gene_name_list[chrom][index1:index2]))
                                        
        mappings_number = self.get_sam_tag_value(mapping.strip().split("\t"),"NH:i:")                         
        #spliting annotation elements by strand
        stranded_genes={}
        for gene_index in genes_on_mapping:
            if gene_list[chrom][gene_index][2] not in stranded_genes:
                stranded_genes[gene_list[chrom][gene_index][2]] = [gene_index]
            else:
                stranded_genes[gene_list[chrom][gene_index][2]].append(gene_index)
                
        #if there are annotation elements in the same strand as mapping then
        #the read is considered only for those annoation elements        
        if strand in stranded_genes:
            for gene_index in stranded_genes[strand]:
                gene_list[chrom][gene_index][5] += \
                                1/int(mappings_number)/len(stranded_genes[strand])
        #if not but there are genes in other strand then the mapping is considered for those
        elif len(stranded_genes) != 0: 
            for gene_strand in stranded_genes:
                for gene_index in stranded_genes[gene_strand]:
                    gene_list[chrom][gene_index][6] += \
                                1/int(mappings_number)/len(stranded_genes[gene_strand])
        #if mappingt does not overlap with any gene       
        else:
            gene_list[chrom][0][7] += 1/int(mappings_number)
        return gene_list
            
    def reads_to_pps4(self,settings,mapping,pp_list,pp_data,pp_pos_list,pp_name_list,
                      start,end,chrom):
        '''
        Looks up pp overlaping with mapping and puts mapping info to gene_list
        '''

        #if pp-s only in one strand but reads in both
        if chrom not in pp_pos_list:
            return pp_data
        
        #get the annotation elements in the positions of the mapping
        ##get positions which are also positions of the annotation elements
        ##and gets directly gene names

        #getting pp-s overlapping with mapping
        index1=bisect.bisect_left(pp_pos_list[chrom],start)
        index2=bisect.bisect(pp_pos_list[chrom],end)
        if index1 == index2:
            pps_on_mapping = pp_name_list[chrom][index1].intersection(pp_name_list[chrom][index1-1])
        else:
            pps_on_mapping = set.union(*map(set,pp_name_list[chrom][index1:index2]))
        #get pp-s which have required overlap with reads
        for pp in pps_on_mapping:
            max_length = max(end-start+1,pp_list[chrom][pp][1]-pp_list[chrom][pp][0])
            overlap = min(end+1,pp_list[chrom][pp][1])-max(start,pp_list[chrom][pp][0])
            pp_name = pp_list[chrom][pp][3]
            #collects total coverage of the pp
            pp_data[pp_name]["total_cov"] += overlap
            #get pp-s which have overlap in suitable range
##            print("A",pp_name)
##            print("B",index1,index2,overlap)
##            print("C",mapping)
##            print("D",start,pp_list[chrom][pp][1],end,pp_list[chrom][pp][1])
            if max_length-settings["non_overlap"] <= overlap <= max_length+settings["non_overlap"]:

                pp_data[pp_name]["count"]["total"][0] += 1 #count
                #collect sequence data
                read_offset = start-pp_list[chrom][pp][0]+settings["non_overlap"]
                #print(start,end,pp_list[chrom][pp][1],pp_list[chrom][pp][0],overlap,max_length,pp_name,read_offset)
                for i, pos in enumerate(range(len(mapping[9]))):
                    pp_data[pp_name]["sequencies"][i+read_offset]["ACGTN".find(mapping[9][i])] += 1
                #uniqness)
                pp_data[pp_name]["uniqness"] += int(self.get_sam_tag_value(mapping,"NH:i:"))
                #rel_cov comes at the end of contig or end of file
                #coverage comes at the end of contig or end of file
                #consensus sequence and consensus quality comes at the end of contig or end of file                
        return pp_data

    def get_sam_tag_value(self,mapping,tag):
        '''
        Gets value of sam tag
        '''
        extra_tags = mapping[11:]
        tag_index = [i for i, j in enumerate(extra_tags) if j.startswith(tag)][0]
        return extra_tags[tag_index].strip(tag)

    def calculate_pp_data(self,settings,pp_data_single):
        '''
        Calculates uniqness, coverage, consensus sequence, consensus quality, genomic sequence
        '''
        #calcualtes only if there are reads
        if pp_data_single["count"]["total"][0] != 0:
            consensus_seq = ""
            consensus_qual = ""
            genomic_seq = ""
            coverage = []
            for i, pos in enumerate(pp_data_single["sequencies"]):
                pos_counts = sum([x for x in pp_data_single["sequencies"][i]])
                coverage.append(str(pos_counts))
                max_value = max(pp_data_single["sequencies"][i])
                #if the most frequent nucleotide is not at least in 50% of the cases give N 
                if pos_counts == 0: #no counts in that position
                    consensus_seq += " "
                    consensus_qual += str(0)
                    genomic_seq += " "
                elif max_value < pos_counts*0.5:
                    consensus_seq += "N"
                    consensus_qual += str(int(round(max_value/pos_counts*10-1,0)))
                    genomic_seq += pp_data_single["raw_genomic_sequence"][i]
                else: # at least in 50% cases
                    #get most frequent nucleotide
                    #(in case of several with equal gives first in the list)
                    #in future give double-degenerate codes as in
                    #http://www.cisred.org/content/methods/help/pfm
                    #print("c")
                    consensus_seq += "ACGTN"[pp_data_single["sequencies"][i].index(max_value)]
                    consensus_qual += str(int(round(max_value/pos_counts*10-1,0)))
                    genomic_seq += pp_data_single["raw_genomic_sequence"][i]
                #adding asterix (*) to the description if some non_overlap allowed
                if settings["non_overlap"] != 0:
                    if i == (settings["non_overlap"]-1) or\
                       i == (len(pp_data_single["sequencies"])-settings["non_overlap"]-1):
                        consensus_seq += "*"
                        consensus_qual += "*"
                        coverage.append("*")
                        genomic_seq += "*"
            #calculate relative coverage
            if settings["non_overlap"] != 0:
                coverage_sum = sum([int(pos) for pos in coverage[settings["non_overlap"]+1:\
                                                            -settings["non_overlap"]-1]])
            else:
                coverage_sum = sum([int(pos) for pos in coverage])      
            pp_data_single["rel_cov"] = round(coverage_sum/pp_data_single["total_cov"],4)
            #calculate uniqness
            pp_data_single["uniqness"]=round(pp_data_single["uniqness"]/\
                                             pp_data_single["count"]["total"][0],2)       
            #reverse and reverse complement all
            if pp_data_single["data"][3] == "-":
                consensus_qual = consensus_qual[::-1]
                coverage = coverage[::-1]
                consensus_seq = self.rev_comp4(consensus_seq)
                genomic_seq = self.rev_comp4(genomic_seq)
            #format coverage
            if pp_data_single["count"]["total"][0] != 0:
                if settings["non_overlap"] == 0:
                    coverage = coverage[0]
                else:
                    coverage = ";".join(coverage[:settings["non_overlap"]*2+1])+\
                               "..."+coverage[settings["non_overlap"]*2+2]+"..."+\
                               ";".join(coverage[-settings["non_overlap"]*2-1:])
            #convert genomic sequence (eg ..*...C.....*A)
            genomic_seq_conv = []
            for j,pos in enumerate(genomic_seq):
                if pos == "*":
                    genomic_seq_conv.append("*")
                elif pos == consensus_seq[i]:
                    genomic_seq_conv.append(".")
                else:
                    genomic_seq_conv.append(pos)
            genomic_seq_conv = "".join(genomic_seq_conv)
            
            pp_data_single["consens_seq"] = consensus_seq
            pp_data_single["consens_qual"] = consensus_qual
            pp_data_single["coverage"], pp_data_single["genomic_seq"] = coverage,genomic_seq_conv

        return pp_data_single

    def rev_comp4(self,seq):
        '''
        Reverse complement sequence
        '''
        seq1 = 'ATCGNRYMKBDHV *TAGCNYRKMVHDB *'
        seq_dict = { seq1[i]:seq1[i+15] for i in range(30) if i < 15}
        return "".join([seq_dict[base] for base in reversed(seq)])

    def write_annotation_statistics(self,gene_list,output_file1,output_file2):
        '''
        Write annotation statistics by libraries
        '''
        f1 = open(output_file1,"w")
        f1.write("############## Annotation by annotation elements ##############\n\n")
        f1.write("Following annotation features has been identified among transcripts:")
        f1.write("STATISTICS BASED ON READ NUMBER:\n\n")
        f1.write("Chromosome\tname\tstart\tend\tstrand\ttype\tsense exon\tsense intron\t"+\
                 "antisense exon\tantisense intron\tundirected features\n")                             
        biotypes = {}
        for chrom in list(gene_list.keys()):
            for gene in gene_list[chrom]:
                #leaves out the number which was gene order number:
                f1.write("\t".join([chrom]+[str(round(x,0)) if isinstance(x, float) else str(x) \
                                            for x in gene[:-1]])+"\n") 
                if gene[4] not in list(biotypes.keys()):
                    biotypes[gene[4]] = gene[5:-1]
                else:
                    biotypes[gene[4]] = [sum(item) for item in zip(gene[5:-1],biotypes[gene[4]])]
                
            #print(gene_list[chrom]) 
        f1.close()
        f2 = open(output_file2,"w")
        f2.write("############## Annotation ##############\n\n")
        f2.write("Following annotation features has been identified among transcripts:")
        f2.write("STATISTICS BASED ON READ NUMBER:\n\n")
        f2.write("biotype\tsense exons\tantisense exons\tundirected features\n")
        for biotype in list(biotypes.keys()):
            f2.write("\t".join([biotype]+[str(round(x,1)) for x in biotypes[biotype]])+"\n")
        f2.close()

    def get_grupped_biotype_reads(self,biotype_counts,exclude_biotypes):
        '''
        Groups values of biotypes which are not listed separately
        '''
        #total value
        groupped_value = {}
        groupped_value = float(biotype_counts["total"])
        for biotype in biotype_counts:
            if biotype == "total":
                continue
            if not biotype in exclude_biotypes:
                continue
            groupped_value -= float(biotype_counts[biotype])
        return groupped_value
                                   
    def write_pp_recalc_file(self,settings,pp_dict,count_type,comment):
        '''
        Writes combined pp spezific quantification or normalization to file.
        '''
        print("\tWrite pp counts:",count_type)
        with open(os.path.join(settings["--output"],"quantify",\
                               "pp_clustered_counts_"+count_type+".tsv"),"w") as f_out:
            #write comment
            f_out.write("#"+comment+"\n")
            #write header
            f_out.write("#pp")
            f_out.write("\tbiotype")
            for library in sorted(settings["libraries"]):
                f_out.write("\t"+library)
            f_out.write("\n")
            #write content
            first_library = sorted(settings["libraries"])[0]
            for pp_name in sorted(pp_dict[first_library]):
                if pp_name == "biotype_counts":
                    continue
                #checks are there multiple values for read count
                #if there was several biotypes then reads can be calculated for both
                for i, x in enumerate(pp_dict[first_library][pp_name]["count"][count_type]):
                    f_out.write(pp_name)
                    if len(pp_dict[first_library][pp_name]["count"][count_type]) == 1:
                        f_out.write("\t"+";".join(pp_dict[first_library][pp_name]["biotype"]))
                    else:
                        f_out.write("\t"+pp_dict[first_library][pp_name]["biotype"][i])
                    for library in sorted(settings["libraries"]):
                        f_out.write("\t"+str(pp_dict[library][pp_name]["count"][count_type][i]))        
                    f_out.write("\n")

    def write_pp_data_file(self,settings,pp_dict,data_type,comment):
        '''
        Writes specific combined pp data to file.
        '''
        print("\tWrite pp data:",data_type)
        with open(os.path.join(settings["--output"],"quantify",\
                               "pp_clustered_"+data_type+".tsv"),"w") as f_out:
            #write comment
            f_out.write("#"+comment+"\n")
            #write header
            f_out.write("#pp")
            f_out.write("\tbiotype")
            for library in sorted(settings["libraries"]):
                f_out.write("\t"+library)
            f_out.write("\n")
            #write content
            first_library = sorted(settings["libraries"])[0]
            for pp_name in sorted(pp_dict[first_library]):
                if pp_name == "biotype_counts":
                    continue
                #checks are there multiple values for read count
                #if there was several biotypes then reads can be calculated for both
                f_out.write(pp_name)
                f_out.write("\t"+";".join(pp_dict[first_library][pp_name]["biotype"]))
                for library in sorted(settings["libraries"]):
                    if pp_dict[library][pp_name]["count"]["total"][0] != 0:
                        f_out.write("\t"+str(pp_dict[library][pp_name][data_type]))
                    else:
                        f_out.write("\t.")
                f_out.write("\n")

    def write_pp_recalc_file_selected(self,settings,pp_dict,count_type,comment):
        '''
        Writes pp-s selected by the minimal count files.
        '''
        print("\tWrite selected pp counts:",count_type)
        first_library = sorted(settings["libraries"])[0]
        max_pp_count = max(max(max(float(y) for y in pp_dict[library][pp_name]["count"][count_type]) for
                               library in settings["libraries"])
                            for pp_name in pp_dict[first_library] if pp_name != "biotype_counts")      
        min_pp_count = min(max(max(float(y) for y in pp_dict[library][pp_name]["count"][count_type]) for
                               library in settings["libraries"])
                            for pp_name in pp_dict[first_library] if pp_name != "biotype_counts")     
        min_range = list(range(len(str(int(float(min_pp_count)))),len(str(int(float(max_pp_count))))))
        f_stat =  open(os.path.join(settings["--output"],"quantify","selected_pps",\
                                   "pp_clustered_stat_"+count_type+".log"),"w")
        f_stat.write("#Number of PP-s with at least given "+count_type+" counts\n")
        f_stat.write("#Count threshold\tpp-s\n")
        f_stat.write("\t".join(["0",str(len(pp_dict[first_library]))])+"\n")
        for minimum in min_range:
            min_count = 10**minimum 
            with open(os.path.join(settings["--output"],"quantify","selected_pps",\
                                   "pp_clustered_counts_"+count_type+"_min_"+str(min_count)+".tsv")
                      ,"w") as f_out:
                #write comment
                f_out.write("#"+comment+". PP-s selected having " +count_type +" counts over "+
                            str(min_count)+"\n")
                #write header
                f_out.write("#pp")
                f_out.write("\tbiotype")
                for library in sorted(settings["libraries"]):
                    f_out.write("\t"+library)
                f_out.write("\n")
                #write content
                pp_count = 0
                for pp_name in sorted(pp_dict[first_library]):
                    if pp_name == "biotype_counts":
                        continue
                    #skips pp-s which are below minimum count
                    pp_max_count = max(max(float(y) for y in pp_dict[library][pp_name]["count"][count_type]) for
                               library in settings["libraries"])
                    if pp_max_count < min_count:
                        continue
                    pp_count += 1
                    #checks are there multiple values for read count
                    #if there was several biotypes then reads can be calculated for both                    
                    for i, x in enumerate(pp_dict[first_library][pp_name]["count"][count_type]):
                        f_out.write(pp_name)
                        if len(pp_dict[first_library][pp_name]["count"][count_type]) == 1:
                            f_out.write("\t"+";".join(pp_dict[first_library][pp_name]["biotype"]))
                        else:
                            f_out.write("\t"+pp_dict[first_library][pp_name]["biotype"][i])
                        for library in sorted(settings["libraries"]):
                            f_out.write("\t"+str(pp_dict[library][pp_name]["count"][count_type][i]))        
                        f_out.write("\n")
            f_stat.write("\t".join([str(min_count),str(pp_count)])+"\n")
        f_stat.close()

    def collect_statistics(self,settings):
        '''
        Collects library specific data to statistics folder
        '''
        tasks = ["trim", "align", "sam_sort", "pseudoSE", "identify"]
        folder = []
        folder.append(os.path.normpath(os.path.join(settings["--input"],"..")))
        folder.append(os.path.normpath(settings["--output"]))
        for library in settings["libraries"]:
            with open(os.path.join(settings["--output"],"quantify","collected_statistics",\
                                   "collected_stat_"+library+".log"),"w") as f_out:
                for task in tasks:
                    for folder_path in folder:
                        if os.path.isdir(os.path.join(folder_path,task)):
                            f_out.write("#"*(len(task)+8)+"\n##  "+task+"  ##\n"+"#"*(len(task)+8)+"\n\n")
                            if task == "identify":
                                for strand in ["For","Rev"]:
                                    f_out.write("#"+strand+"\n")
                                    with open(os.path.join(folder_path,task,task+"_info",
                                                           library+"_"+strand+"_"+task+"info.log")) as f_in:
                                        for line in f_in:
                                            f_out.write(line)
                                        f_out.write("\n")                                
                            else:
                                with open(os.path.join(folder_path,task,task+"_info",
                                                       library+"_"+task+"info.log")) as f_in:
                                    for line in f_in:
                                        f_out.write(line)
                                    f_out.write("\n")
                            break
                        
    def test_errors(self,settings):
        '''
        Testing size of output files
        '''
        for data_type in ["consens_seq","genomic_seq","consens_qual","coverage","uniqness","rel_cov"]:
            if not os.path.isfile(os.path.join(settings["--output"],"quantify",\
                                            "pp_clustered_"+data_type+".tsv")):
                sys.exit('Task "quantify" incomplete!A')
            if os.stat(os.path.join(settings["--output"],"quantify",\
                                   "pp_clustered_"+data_type+".tsv")).st_size == 0:
                sys.exit('Task "quantify" incomplete!B')
        for count_type in ["total","RPM","biotype_RPM","groupped_biotype_RPM"]:
            if not os.path.isfile(os.path.join(settings["--output"],"quantify",\
                                        "pp_clustered_counts_"+count_type+".tsv")):
                sys.exit('Task "quantify" incomplete!C')		
            if os.stat(os.path.join(settings["--output"],"quantify",\
                                   "pp_clustered_counts_"+count_type+".tsv")).st_size == 0:
                sys.exit('Task "quantify" incomplete!D')
