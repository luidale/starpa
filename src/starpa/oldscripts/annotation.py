#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Converts GFF and GFF3 format annotation to BED annotation form.
Groups different ncRNA-s as ncRNAs.
Excludes various annotation elements.
Skips annotation region description line.
'''

import copy

class annotation():
    def __init__(self,settings):
        self.settings = settings
        self.settings = self.convert(settings)

    def return_settings(self):
        return self.settings
        
    def convert(self,settings):
        '''
        Script doing the annoation
        '''
        print("Convert annotation file")
        annotations = {} #dictionary for annotations
        biotypes=set() #set with identified biotypes
        #open input_gff_file
        f_in = open(settings["quantify"]["quantify_annotation_file"])
        line = f_in.readline()

        for next_line in f_in:
            #skip headers and comments and first line of the region,
            #which cover the whole region
            while line[0] == "#":
                line = copy.deepcopy(next_line)
                next_line = f_in.readline()
                
            #skip non-interesting annotations
            line,next_line = self.skip_lines(settings,line,next_line,f_in)

            #Get gbkey
            gb_key = self.get_annotation_tag(line.strip().split("\t")[8],"gbkey")
            
            #get single line interesting annotations
            if gb_key in settings["quantify"]["quantify_single_line_elements"]:
                annotations = self.add_one_lane_element(line.strip().split("\t"),annotations,gb_key)
                line = copy.deepcopy(next_line)
                
            #genes
            elif gb_key == "Gene":
                #get line group of the gene
                group,next_line = self.get_line_group(settings,line,next_line,f_in)
                
                ##get data
                #main data is extracted from first line of the group which is gene
                group_line = group[0]
                group_line_split = group_line.strip().split("\t")
                chrom = group_line_split[0]
                biotype = self.get_annotation_tag(group_line_split[8],"gene_biotype")
                biotypes.add(biotype)
                name = self.get_annotation_tag(group_line_split[8],"Name")
                description = self.get_annotation_tag(group_line_split[8],"description")
                product = self.get_annotation_tag(group_line_split[8],"product")
                note = self.get_annotation_tag(group_line_split[8],"Note")
                
                #collect exons
                exons = set()
                #parsing lines after the first gene line
                for group_line in group[1:]:
                    #skips line "mRNA"
                    if group_line.strip().split("\t")[2] == "mRNA":
                        continue

                    #gets description and product if still missing
                    if not description:
                        description = self.get_annotation_tag(group_line.strip().split("\t")[8]\
                                                         ,"description")
                    if not product:
                        product = self.get_annotation_tag(group_line.strip().split("\t")[8],\
                                                          "product")
                    if not note:
                        note = self.get_annotation_tag(group_line.strip().split("\t")[8],"Note")

                    #collect annotation types    
                    group_line_split = group_line.strip().split("\t")
                    if biotype in settings["quantify"]["quantify_ncRNA_types"]:
                        exons.add((int(group_line_split[3]),int(group_line_split[4]),\
                                   group_line_split[6],name,"ncRNA",biotype,\
                                   str(description),str(product),str(note)))
                    else:
                        exons.add((int(group_line_split[3]),int(group_line_split[4]),\
                                   group_line_split[6],name,biotype,\
                                   str(description),str(product),str(note)))
                
                if chrom not in annotations:
                    annotations[chrom] = []
                    
                #no exons given
                if len(exons) == 0:
                    """In E.coli some pseudogenes are given as two separate
                    genes when it is distrupted by prophage or something
                    similar. First part of the pseudogene contains only gene
                    description and no CDS while the second gene contains both
                    parts as CDS.
                    To avoid repetition of CDS in BED file it is checked does
                    the gene without exons has the same name as the following
                    line. If following line has the same name annotation tag
                    then the first is skipped. Other genes without exons
                    (single line annotations) are still added to annotations"""
                    if len(next_line.strip().split("\t")) > 8:
                        if self.get_annotation_tag(line.strip().split("\t")[8],"Name") == \
                            self.get_annotation_tag(next_line.strip().split("\t")[8],"Name"):
                            continue
                        else: #single line annotation
                            annotations = self.add_one_lane_element(\
                                            line.strip().split("\t"),annotations,biotype)
                        
                else: #exon(s) given
                    exons = sorted(list(exons))
                    exons[0]= list(exons[0])
                    for i, exon in enumerate(exons):
                        annotations[chrom].append(list(exon))

##                      NOT IN USE
##                        #If splicing is allowed then introns are looked for
##                        if splicing:
##                            #testing is there introns
##                            if i != len(exons)-1: #testing is it last exon
##                                #testing gap between exons
##                                if int(exons[i+1][0])-int(exons[i][1]) > 0: 
##                                    annotations[chrom].append([exons[i][1],exons[i+1][0],\
##                                                               exons[i][2],name,"intron",\
##                                                               exons[i][4],exons[i][5]])

            #prints unrecognized annotation elements
            else:
                print("Unrecognized annotation elements:",gb_key)
                
            line = copy.deepcopy(next_line)
        f_in.close()


        #write gene biotypes and the annotation types which were considered as ncRNA
        ncRNA_types_found = set(biotypes).intersection(set(\
                                settings["quantify"]["quantify_ncRNA_types"]))
        other_biotypes = set(biotypes).difference(ncRNA_types_found)
        annotation_convert_info_file = settings["quantify"]["quantify_annotation_file"]+".conv_info"
        with open(annotation_convert_info_file,"w") as f_info:
            f_info.write("gbkeys:\n")
            f_info.write("\tSkipped:"+" , ".join(sorted(settings["quantify"]\
                                                        ["quantify_keys_to_skip"]))+"\n")
            f_info.write("\tIncluded:"+" , ".join(sorted(settings["quantify"]\
                                                         ["quantify_single_line_elements"]+\
                                                         ["Gene"]))+"\n")
            f_info.write("\t\tGene biotypes in output:"+\
                                             " , ".join(sorted(list(other_biotypes)+["ncRNA"]))+"\n")
            f_info.write("\t\tGene biotypes groupped as ncRNA:"+\
                                             " , ".join(sorted(list(ncRNA_types_found)))+"\n")


        #write output
        self.write_annotation(settings,annotations)
            
    def skip_lines(self,settings,line,next_line,f_in):
        '''
        Skips lines of input files with annotation elements which are meant to be skipped.
        '''
        #remove non-interesting annotations
        while self.get_annotation_tag(line.strip().split("\t")[8],"gbkey") in \
                                      settings["quantify"]["quantify_keys_to_skip"]:
            line = copy.deepcopy(next_line)
            next_line = f_in.readline()
            
            #remove header and comment lines
            while line[0] == "#":
                line = copy.deepcopy(next_line)
                next_line = f_in.readline()       
            #testin last line
            if line.strip() =="":
                break
        return line, next_line

    def get_annotation_tag(self,string,tag):
        '''
        Gets annotation tag
        '''
        #print(string,tag)
        for element in string.split(";"):
            if element.startswith(tag+"="):
                return element[len(tag)+1:]
        return False

    def add_one_lane_element(self,line_split,annotations,annotation_type):
        '''
        Adds elements which consists form single line
        '''
        #check chromosome/contig
        if line_split[0] not in annotations:
            annotations[line_split[0]]=[]
            
        name = self.get_annotation_tag(line_split[8],"Name")
        description = self.get_annotation_tag(line_split[8],"Note")
        product = self.get_annotation_tag(line_split[8],"product")
        if annotation_type == "repeat_region":
            if not name:
                try:
                    name = self.get_annotation_tag(line_split[8],"Note").split(" ")[0]
                except AttributeError:
                    #exception for trypanosoma repeat_regions
                    name = self.get_annotation_tag(line_split[8],"ID")
        elif annotation_type == "telomere":
            if not name:
                name = self.get_annotation_tag(line_split[8],"Note").split("%")[0]
        elif annotation_type == "rep_origin":
            if not name:
                name = self.get_annotation_tag(line_split[8],"Note").split("~")[0]
                if len(self.get_annotation_tag(line_split[8],"Note").split("~")) == 1:
                    name = self.get_annotation_tag(line_split[8],"Note").split("%")[0]
        annotations[line_split[0]].append([int(line_split[3]),int(line_split[4]),\
                                           line_split[6],str(name),annotation_type,\
                                           str(description),str(product)])
        return annotations
    
    def get_line_group(self,settings,line,next_line,f_in):
        '''
        Gets line groups which belong to the same gene.
        Stops when file ends, comment line, new gene
        '''
        group = [line]
        while True:
            if next_line.strip() == "":
                break
            elif next_line[0] == "#":
                break
            elif self.get_annotation_tag(next_line.strip().split("\t")[8],"gbkey") in \
                                            settings["quantify"]["quantify_keys_to_skip"]:
                break
            elif self.get_annotation_tag(next_line.strip().split("\t")[8],"gbkey") in \
                                            settings["quantify"]["quantify_single_line_elements"]:
                break
            elif self.get_annotation_tag(next_line.strip().split("\t")[8],"gbkey") == "Gene":
                break

            group.append(next_line)
            line = copy.deepcopy(next_line)
            next_line = f_in.readline()
        return group, next_line
        
    def write_annotation(self,settings,annotations):
        '''
        Writes annotation file in BED format
        '''
        annotation_BED_file = ".".join(settings["quantify"]["quantify_annotation_file"].\
                                       split(".")[:-1])+".BED"
        f_out_BED = open(annotation_BED_file,"w")
        #f_out_SAF = open(out_file+name_sufix+".SAF","w")
        for chrom in sorted(annotations):
            for annotation in sorted(annotations[chrom]):
                #writes on one line:
                #chrom, start,end,strand,name, gene_biotype, biotype;product;description
                f_out_BED.write("\t".join([chrom]+[str(annotation[0]-1),str(annotation[1]),\
                                        annotation[3],"1",annotation[2],\
                                        annotation[4]]+[";".join(sorted(\
                                            list(set([x for x in annotation[5:] if x!="False"])),\
                                                key=len))])+"\n")
                #f_out_SAF.write("\t".join([annotation[3]]+[chrom]+\
                #[str(x) for x in annotation[:3]])+"\n")
        f_out_BED.close()
        #f_out_SAF.close()
        settings["quantify"]["quantify_annotation_file"] = annotation_BED_file
