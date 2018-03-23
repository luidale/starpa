#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Usage:
    starpa [-hv]
    starpa -s <start_task> -e <end_task> -c <config_file> -i <input> -o <output>
    
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
  
'''


from docopt import docopt
from schema import Schema, And, Or, Use, SchemaError
from starpa.__init__ import __version__

from starpa.trim import trim
from starpa.align import align
from starpa.sort import sam_sort
from starpa.pseudoSE import pseudoSE
from starpa.identify import identify
from starpa.cluster import cluster
from starpa.quantify import quantify

import sys
import os
import math
from shutil import copy2

tasks = ["trim", "align", "sam_sort", "pseudoSE", "identify", "cluster", "quantify"]

def get_tasks(args,tasks):
    '''
    Gets tasks to be run
    '''
    tasks_to_run = []
    if "--start" in args:
        tasks_to_run = tasks[tasks.index(args["--start"]):]
    else:
        tasks_to_run = tasks
    if "--end" in args:
        tasks_to_run = tasks_to_run[:tasks_to_run.index(args["--end"])+1]
    return tasks_to_run

def check_arguments(args, tasks):
    #test config file and tasks
    schema = Schema({
        "--start": And(lambda s: s in set(tasks), error="Invalid start task"),
        "--end": And(lambda s: s in set(tasks), error="Invalid end task"),
        "--config": And(open, error='<config_file> should be readable'),
        "--input": And(os.path.isdir, error='<input_folder> should be folder'),
        object: object })
    try:
        args = schema.validate(args)
    except SchemaError as e:
        exit(e)

    #test of start and end tasks match   
    if tasks.index(args["--start"]) > tasks.index(args["--end"]):
        sys.exit("Check tasks. End task earlier than start task")

    #read in config file
    args = read_in_arguments(args,tasks)

    #test general arguments
    schema = Schema({
        "library_file": And(open, error="'library_file' should be readable"),
        "genome": And(open, error="'genome' should be readable"),
        "paired": And(lambda s: s in {"True","False"},(Use(str_to_bool)),\
                                   error="'paired' should be Boolean (True or False)"),
        "samtools_call": And(Use(str), error="'samtools_call' should be string"),
        "bedtools_call": And(Use(str), error="'bedtools_call' should be string"),
##        "featureCounts_call": And(Use(str), error="'featureCounts_call' should be string"),
        "cd_hit_est_call": And(Use(str), error="'cd_hit_est_call' should be string"),
        "CPUs": And(Use(int), (lambda s: s > 0), error="'CPUs' should be integer"),
        "min_length": And(Use(int),(lambda s: s > 0),\
                                    error="'min_length' should be positive integer"),
        "max_length": And(Use(int),(lambda s: s > 0),\
                                    error="'max_length' should be positive integer"),
##        "overlap_fraction": And(Use(float),(lambda s: 1 >= s >= 0),\
##                                    error="'max_length' should be positive integer"),
        "non_overlap": And(Use(int),(lambda s: s > 0),\
                                    error="'non_overlap' should be positive integer"),
        "min_pp_reads": And(Use(int),(lambda s: s > 0),\
                                    error="'min_pp_reads' should be positive integer"),
		"samtools_threads": And(Use(int),(lambda s: s > 0),\
                                    error="'samtools_thread' should be positive integer"),
        object: object
        })
    try:
        args = schema.validate(args)
    except SchemaError as e:
        exit(e)
    #compare min_length and max_length
    if args["min_length"] > args["max_length"]:
        sys.exit("'min_length' should not be bigger than 'max_length'")
        
    #get tasks
    tasks_to_run = get_tasks(args,tasks)

    #test arguments by tasks
    for task in tasks_to_run:
        if task == "trim":
            schema = Schema({
                "trim_call": And(str, error="trim_call should be string"),
                "trim_min_qual": And(Use(int), error="'trim_quality' should' be integer"),
                "trim_overlap": And(Use(int), error="'trim_overlap' should' be integer"),
                "trim_adapter_for": And(open, error="'trim_adapter_for' should be readable"),
                "trim_adapter_rev": And(open, error="trim_adapter_rev' should be readable"),
                "trim_threads": And(Use(int),(lambda s: s == 1), \
                                   error="'trim_threads' should be integer 1"),
                "trim_quality_base": And(lambda s: s in {"33","64"},\
                                         error="'trim_quality_base should be 33 or 64"),
                "trim_input_file_suffix_for": And(str, \
                                        error="'trim_input_file_suffix_for' should be string"),
                "trim_input_file_suffix_rev": And(str, \
                                        error="'trim_input_file_suffix_rev' should be string"),
                "trim_input_file_suffix_SE": And(str, \
                                        error="'trim_input_extention_SE' should be string")
                })

        elif task in {"align","sam_sort"}:
            schema = Schema({        
                "align_call": And(str, error="'align_call' should be string"),
                "align_index_call": And(str, error="'align_index_call' should be string"),
                "align_threads": And(Use(int),(lambda s: s > 0),\
                                     error="'align_threads' should be positive integer"),
                "align_quality_base": And(lambda s: s in {"33","64"},\
                                    error="'align_quality_base' should be 33 or 64"),
                "align_sensitive": And(lambda s: s in {"True","False",True,False},\
                                       (Use(str_to_bool)),\
                                    error="'align_sensitive' should be Boolean (True or False)"),
                "align_input_file_suffix_for": And(str, \
                                    error="'align_input_file_suffix_for' should be string"),
                "align_input_file_suffix_rev": And(str, \
                                    error="'align_input_file_suffix_rev' should be string"),
                "align_input_file_suffix_SE": And(str, \
                                    error="'align_input_extention_SE' should be string")
                })

        elif task == "pseudoSE":
            schema = Schema({        
                "pseudoSE_input_file_suffix": And(str,\
                        error="'pseudoSE_input_file_suffix' should be string"),
                "pseudoSE_max_mappings": And(Use(int),\
                        error="'pseudoSE_max_mappings' should be integer"),
                "pseudoSE_oligoA": And(lambda s: s in {"True","False",True,False},\
                                       (Use(str_to_bool)),\
                        error="'pseudoSE_oligoA' should be Boolean (True or False)"),
                "pseudoSE_quality_base": And(lambda s: s in {"33","64"},\
                                    error="'pseudoSE_quality_base' should be 33 or 64"),
                "pseudoSE_allowed_mismatch": And(Use(int),\
                        error="'pseudoSE_allowed_mismatch' should be integer"),
                "pseudoSE_mismatch_precentage": And(Use(int),(lambda s: s in range(0,100+1)),\
                        error="'pseudoSE_mismatch_precentage' should be integer in range 0-100"),
##               "pseudoSE_max_read_length": And(Use(int),\
##                        error="'pseudoSE_max_read_length' should be integer")
##                "pseudoSE_sorted_by_pos": And(lambda s: s in {"True","False",True,False},\
##                                       (Use(str_to_bool)),\
##                        error="'pseudoSE_sorted_by_pos' should be Boolean (True or False)")
                })
            
##        if task == "sam_sort":
##            schema = Schema({
##                "sam_sort_sorted_by_pos": And(lambda s: s in {"True","False",True,False},\
##                                       (Use(str_to_bool)),\
##                        error="'pseudoSE_sorted_by_pos' should be Boolean (True or False)"),
##                object: object
##                })

        elif task == "identify":
            schema = Schema({   
                "identify_flaimapper_parameters": And(open, \
                        error="'identify_flaimapper_parameter' should be readable"),
##                "identify_overlap": And(Use(int),\
##                        error="'identify_overlap' should be integer"),
                "identify_flaimapper_call": And(Use(str),\
                        error="'identify_flaimapper_call' should be string"),
                "identify_input_file_suffix": And(str,\
                        error="'identify_input_file_suffix' should be string"),
                "identify_split_step": And(Use(int),\
                        error="'identify_split_step' should be integer")
                })

        elif task == "cluster":
            schema = Schema({   
                "cluster_min_contig_length": And(Use(int),(lambda s: s > 0), \
                        error="'cluster_min_contig_length' sshould be positive integer"),
                "cluster_min_contig_cov": And(Use(int),(lambda s: s > 0), \
                        error="'cluster_min_contig_cov' should be positive integer"),
                "cluster_min_contig_reads": And(Use(int),(lambda s: s > 0), \
                        error="'cluster_min_contig_reads' should be positive integer"),
                "cluster_min_contig_length_meta": And(Use(int),(lambda s: s > 0), \
                        error="'cluster_min_contig_length_meta' should be positive integer"),
                "cluster_min_contig_cov_meta": And(Use(int),(lambda s: s > 0), \
                        error="'cluster_min_contig_cov_meta' should be positive integer"),
                "cluster_min_contig_reads_meta": And(Use(int),(lambda s: s > 0), \
                        error="'cluster_min_contig_reads_meta' should be positive integer"),
                "cluster_rel_cov_list": And(Use(str_to_list_int),\
                        [int,float], error="'cluster_rel_cov_list' should be a list as "+\
                                            "[X] or [X,...,Y] containing numbers below 1"),
                "cluster_rel_cov_size_range": And(Use(str_to_list_int),\
                        [int,float], error="'cluster_rel_cov_size_range' should be a list as "+\
                                            "[X] or [X,...,Y] containing integers bigger than 0"),
                "cluster_input_file_suffix": And(str,\
                        error="'cluster_input_file_suffix' should be string"),
                object: object
                })

        elif task == "quantify":
            schema = Schema({
##                "quantify_bam_file_location": And(os.path.isdir,\
##                        error="'quantify_bam_file_location' should be folder"),
##                "quantify_sam_file_location": And(os.path.isdir,\
##                        error="'quantify_sam_file_location' should be folder"),
                "quantify_sam_file_suffix": And(str,\
                        error="'quantify_sam_file_suffix' should be string"),
                "quantify_annotation_file": And(open,\
                        error="'quantify_annotation_file' should be readable"),
                "quantify_keys_to_skip": And(Use(str_to_list_str),\
                        [str], error="'quantify_keys_to_skip' should be a list as "+\
                                            "[X] or [X,...,Y] containing strings"),
                "quantify_single_line_elements": And(Use(str_to_list_str),\
                        [str], error="'quantify_single_line_elements' should be a list as "+\
                                            "[X] or [X,...,Y] containing strings"),
                "quantify_ncRNA_types": And(Use(str_to_list_str),\
                        [str], error="'quantify_ncRNA_types' should be a list as "+\
                                            "[X] or [X,...,Y] containing strings"),
                "quantify_non_groupped_biotypes": And(Use(str_to_list_str),\
                        [str], error="'quantify_non_groupped_biotypes' should be a list as "+\
                                            "[X] or [X,...,Y] containing strings"),
                object: object
                })
            
        try:
            if task in args:
                args[task] = schema.validate(args[task])
            if task == "sam_sort":
                args["align"] = schema.validate(args["align"])
        except SchemaError as e:
            exit(e)

        if task in {"identify","quantify"}:
        #adds read-pp overlap ans size range parametes
            args["size_range"],args["overlap_range"] =\
                                            get_range(args["min_length"],\
                                                      args["max_length"],\
                                                      args["non_overlap"])
        if task == "cluster":
            if len(args["cluster"]["cluster_rel_cov_list"]) !=\
               len(args["cluster"]["cluster_rel_cov_size_range"])+1:
                sys.exit('"cluster_rel_cov_list" should be longer by one element\
                         than "cluster_rel_cov_range')

        #test that sma_sort is nat start task when alignment is in
        #sensitive mode
        if args["--start"] == "sam_sort" and args["align"]["align_sensitive"]:
            sys.exit('Check tasks or "align_sensitive" parameter.\n'+\
                     '"sam_sort" can not be start task while align is in '+\
                     'sensitive mode: "align_sensitive" = True')
        
    return args,tasks_to_run    

def get_range(min_length,max_length,non_overlap):
    '''
    Calculates length and overlap ranges for bedtools intersect.
    Overlap value for given range will fullfill non-overlap for
    given number of nucleotides
    '''
    size_range = []
    overlap = []
    #if overlap between pp and read has to be exact
    if non_overlap == 0:
        overlap.append(non_overlap)
        return size_range,overlap

    #gets first overlap
    new_overlap = rounddown((min_length-non_overlap)/(min_length))
    ##if overlap value does not result wanted non-overlap it will be
    ##increased stepwise and tested again
    while not (min_length*new_overlap <= max(min_length-non_overlap,min_length)) \
          or not (min_length/new_overlap <= min_length+non_overlap+1):
        new_overlap = new_overlap*1.002
        if new_overlap > 1:
            sys.exit("Decrease 'non_overlap' parameter. It can be\
                     upto around square root of the parameter 'min_length' \
                     (the minimal length of the reads)")
    
    overlap.append(new_overlap)

    #scan whole read length range
    for x in range(min_length+1,max_length+1):
        #overlap parameter suits with the current length
        if (x*overlap[-1] < x-non_overlap) and (x/overlap[-1] < x+non_overlap+1):
            continue

        #overlap parameter does not suit with the current length
        else:     
            new_overlap = rounddown((x-non_overlap)/(x))
            ##if overlap value does not result wanted non-overlap it will be
            ##increased stepwise and tested again
            while not (x*new_overlap <= max(x-non_overlap,min_length)) \
                  or not (x/new_overlap <= x+non_overlap+1):
                new_overlap = new_overlap*1.002
                if new_overlap > 1:
                    sys.exit("Decrease 'non_overlap' parameter. It can be\
                         upto around square root of the parameter 'min_length' \
                         (the minimal length of the reads)")
                    
            overlap.append(new_overlap)
            size_range.append(x)
    return size_range,overlap

def rounddown(x):
    '''
    Rounds the value down
    '''
    return int(math.floor(x * 100000.0)) / 100000.0

def read_in_arguments(args,tasks):
    with open(args["--config"]) as f_in:
        for line in f_in:
            if line.strip().startswith("#"):
                continue
            if line.strip().strip("\t") == "":
                continue
            #get arguments string
            if line.find("#") != -1:
                line = line.strip().split("#")[0].strip("\t")
            else:
                line = line.strip()
            #get arguments
            ## by tasks
            if line.split("_")[0] in tasks:
                task = line.strip().split("_")[0]
                argument,argument_value = \
                                (x for x in line.strip().split(" = "))
                if task not in args:
                    args[task] = {}
                args[task][argument] = argument_value
            ##without tasks
            else:
                args[line.split(" = ")[0]] = line.split(" = ")[1]
    return args

def check_input_files(args):
    '''
    Checks existence of suitable files in input folder.
    File have to be bigger than 0.
    '''

    strand_list = ["For", "Rev"]
    
    if args["--start"] == "trim":
        for file in os.listdir(args["--input"]):
            if args["paired"]:        
                if file.endswith(args["trim"]["trim_input_file_suffix_for"]):
                    if os.stat(os.path.join(args["--input"],file)).st_size != 0:
                        break
            else:
                if file.endswith(args["trim"]["trim_input_file_suffix_SE"]):
                    if os.stat(os.path.join(args["--input"],file)).st_size != 0:
                        break               
        else:
            if args["paired"]:
                sys.exit('Suitable files (ending with "'+\
                         args["trim"]["trim_input_file_suffix_for"]+\
                         '") missing in input folder')
            else:
                sys.exit('Suitable files (ending with "'+\
                         args["trim"]["trim_input_file_suffix_SE"]+\
                         '") missing in input folder')
                                
    elif args["--start"] == "align":
        for file in os.listdir(args["--input"]):
            if args["paired"]:        
                if file.endswith(args["align"]["align_input_file_suffix_for"]):
                    if os.stat(os.path.join(args["--input"],file)).st_size != 0:
                        break
            else:
                if file.endswith(args["align"]["align_input_file_suffix_SE"]):
                    if os.stat(os.path.join(args["--input"],file)).st_size != 0:
                        break               
        else:
            if args["paired"]:
                sys.exit('Suitable files (ending with "'+\
                         args["align"]["align_input_file_suffix_for"]+\
                         '") missing in input folder')
            else:
                sys.exit('Suitable files (ending with "'+\
                         args["align"]["align_input_file_suffix_SE"]+\
                         '") missing in input folder')
            
    elif args["--start"] == "sam_sort":
        for file in os.listdir(args["--input"]):
            if file.endswith(".sam"):
                if os.stat(os.path.join(args["--input"],file)).st_size != 0:
                    break
        else:
            sys.exit('Suitable files (ending with ".sam") missing in input folder')

    elif args["--start"] == "pseudoSE":
        for file in os.listdir(args["--input"]):
            if file.endswith(args["pseudoSE"]["pseudoSE_input_file_suffix"]):
                if os.stat(os.path.join(args["--input"],file)).st_size != 0:
                    break
        else:
            sys.exit('Suitable files (ending with "'+\
                         args["pseudoSE"]["pseudoSE_input_file_suffix_SE"]+\
                         '") missing in input folder')

    elif args["--start"] == "identify":
        for file in os.listdir(args["--input"]):
            if file.endswith(args["identify"]["identify_input_file_suffix"]):
                if os.stat(os.path.join(args["--input"],file)).st_size != 0:
                    break
        else:
            sys.exit('Suitable files (ending with "'+\
                         args["identify"]["identify_input_file_suffix"]+\
                         '") missing in input folder')
            
    elif args["--start"] == "cluster":
        #test existing of bam folder
        if not os.path.isdir(os.path.join(args["--input"],"bam")):
            sys.exit("<input_folder> should contain folder 'bam'")
        #test existence of proper bam files
        bam_files = os.listdir(os.path.join(args["--input"],"bam"))
        for library in args["libraries"]:
            for strand_name in strand_list:
                if library+"_"+strand_name+".bam" not in bam_files:
                    sys.exit("File "+library+"_"+strand_name+".bam missing from "+
                             "folder "+os.path.join(args["--input"],"bam"))
        #test existence of proper counted pp files
        bed_files = os.listdir(os.path.join(args["--input"]))
        for library in args["libraries"]:
            for strand_name in strand_list:
                if library+"_"+strand_name+\
                        args["cluster"]["cluster_input_file_suffix"] not in bed_files:
                    sys.exit("File "+library+"_"+strand_name+\
                        args["cluster"]["cluster_input_file_suffix"]+" missing from "+
                             "folder "+os.path.join(args["--input"]))

        #test existence of folder for pseudoSE info
        if not os.path.isdir(os.path.normpath(os.path.join(args["--input"],\
                                                args["cluster"]["cluster_pseudoSE_location"]))):
                        sys.exit(os.path.normpath(os.path.join(args["--input"],\
                                    args["cluster"]["cluster_pseudoSE_location"]))+\
                     ' should be a folder')
        #test existence of pseudoSE info files
        info_files = os.listdir(os.path.normpath(os.path.join(args["--input"],\
                                    args["cluster"]["cluster_pseudoSE_location"])))
        for library in args["libraries"]:
            if library+"_pseudoSEinfo.log" not in info_files:
                sys.exit("File "+library+"_pseudoSE.log missing from "+
                         "folder "+os.path.normpath(os.path.join(args["--input"],\
                                    args["cluster"]["cluster_pseudoSE_location"])))

    elif args["--start"] == "quantify":
        #testing clustered pp file
        if not os.path.isfile(os.path.join(args["--input"],"pp_metacontig.BED")):
            sys.exit('File "pp_metacontig.BED" missing in input folder '+
                     os.path.join(args["--input"]))
##        #test existence of bam folder
##        if not os.path.isdir(os.path.normpath(os.path.join(args["--input"],\
##                                    args["quantify"]["quantify_bam_file_location"]))):
##            sys.exit(os.path.normpath(os.path.join(args["--input"],\
##                                    args["quantify"]["quantify_bam_file_location"]))+\
##                     ' should be a folder')                  
##        #test existence of proper bam files
##        bam_files = os.listdir(os.path.normpath(os.path.join(args["--input"],\
##                                    args["quantify"]["quantify_bam_file_location"])))
##        for library in args["libraries"]:
##            for strand_name in strand_list:
##                if library+"_"+strand_name+".bam" not in bam_files:
##                    sys.exit("File "+library+"_"+strand_name+".bam missing from "+
##                             "folder "+os.path.normpath(os.path.join(args["--input"],\
##                                    args["quantify"]["quantify_bam_file_location"])))
        #test existence of sam folder
        if not os.path.isdir(os.path.normpath(os.path.join(args["--input"],\
                                    args["quantify"]["quantify_sam_file_location"]))):
            sys.exit(os.path.normpath(os.path.join(args["--input"],\
                                    args["quantify"]["quantify_sam_file_location"]))+\
                     ' should be a folder')
        #test existence of proper sam files
        sam_files = os.listdir(os.path.normpath(os.path.join(args["--input"],\
                                    args["quantify"]["quantify_sam_file_location"])))
        for library in args["libraries"]:
            if library+args["quantify"]["quantify_sam_file_suffix"] not in sam_files:
                sys.exit("File "+library+args["quantify"]["quantify_sam_file_suffix"]+\
                         " missing from folder "+os.path.normpath(os.path.join(args["--input"],\
                                    args["quantify"]["quantify_sam_file_location"])))
            
                   

def check_output_folder(args):
    '''
    Creates and checks output folder
    '''
    #make output folder if needed
    if not os.path.exists(args["--output"]):
        os.makedirs(args["--output"])
    if not os.path.exists(os.path.join(args["--output"],"parameters")):
        os.makedirs(os.path.join(args["--output"],"parameters"))
    #test output parameter
    schema = Schema({
        "--output": And(os.path.isdir, error='<output_folder> should be folder'),
        object: object })
    try:
        args = schema.validate(args)
    except SchemaError as e:
        exit(e)
    return args
         
def get_libraries(args):
    '''
    Reads library names into arguments
    '''
    #create list of libraries
    with open(args["library_file"]) as f:
        #list comprehenisons, removes header
        libraries_dic = {line.strip().split("\t")[0]: line.strip().split("\t")[1:] for line\
                     in f if line[0] != "#"}
        args["libraries"] = libraries_dic
    return args

def str_to_bool(s):
    '''
    Converts strings "True" and "False" to bool
    '''
    if s == 'True':
         return True
    elif s == 'False':
         return False
    else:
        return s

def str_to_list_int(s):
    '''
    Converts list given as a string into real list with integers
    '''

    if (s[0] != "[") or (s[-1] != "]"):
        return s
    s = s.strip("[").strip("]").split(",")
    s_list = []
    for x in s:
        try:
            if len(x.split(".")) > 1:
                s_list.append(round(float(x),len(x.split(".")[1])))  
            else:
                s_list.append(int(x))
        except:
            return s
            
    return s_list

def str_to_list_str(s):
    '''
    Converts list given as a string into real list with integers
    '''

    if (s[0] != "[") or (s[-1] != "]"):
        return s
    s = s.strip("[").strip("]").split(",")
    #strip spaces out
    s = [x.strip() for x in s]    
    return s    

def write_arguments(args):
    '''
    Writes arguments in to file in output folder.
    Also copies argument file into output folder.
    '''
    #write arguments
    with open(os.path.join(args["--output"],"parameters","arguments.txt"),"w") as f_out:
        for argument in args:
            if argument not in ["--help","--version"]:
                if argument.startswith("--"):
                    f_out.write("\t".join([argument,str(args[argument])])+"\n")

    #copy config file
    copy2(args["--config"],os.path.join(args["--output"],"parameters"))
    #copy library file
    copy2(args["library_file"],os.path.join(args["--output"],"parameters"))
    #copy adapterfiles
    copy2(args["identify"]["identify_flaimapper_parameters"],os.path.join(args["--output"],"parameters"))


########
#SCRIPT#
########

def main():    
    args = docopt(__doc__,version=__version__)
    args,tasks_to_run = check_arguments(args, tasks)
    args = get_libraries(args)
    check_input_files(args)
    args = check_output_folder(args)
    write_arguments(args)



    #Run tasks
    if "trim" in tasks_to_run:
        print("Trim")
        trim(args)

    if "align" in tasks_to_run:
        #sensitive mode includes task sort
        if args["align"]["align_sensitive"]:
            print("Align I, L 22")
            align(args,"I",tasks_to_run[0])           
            print("Sort SAM I")
            sam_sort(args,"I",tasks_to_run[0])
            print("Align II, L 14")
            align(args,"II",tasks_to_run[0])
            print("Sort SAM II")
            sam_sort(args,"II",tasks_to_run[0])
        else:
            print("Align, L 22")
            align(args,"all",tasks_to_run[0])           
                  
    ##If align in sensitive mode then sam_sort is already done
    if (not args["align"]["align_sensitive"]) or (tasks_to_run[0] == "sam_sort"):
        if "sam_sort" in tasks_to_run:
            print("Sort SAM")
            sam_sort(args,"all",tasks_to_run[0])
            
    if "pseudoSE" in tasks_to_run:
        print("PseudoSE")
        pseudoSE(args,tasks_to_run[0])

    if "identify" in tasks_to_run:
        print("Identify pp-s")
        identify(args,tasks_to_run[0])

    if "cluster" in tasks_to_run:
        print("Cluster pp-s")
        cluster(args,tasks_to_run[0])
            
    if "quantify" in tasks_to_run:
        print("Quantify pp-s")
        quantify(args,tasks_to_run[0])
                             
if __name__ == '__main__':
    main()
