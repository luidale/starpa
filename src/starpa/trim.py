#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
To run read trimming with cutadapt.

To do:

'''

import os
import sys
import multiprocessing as mp

class trim():
    def __init__(self,settings):
        self.settings = settings
        self.make_folder(settings)
        if settings["paired"]:
            pool = mp.Pool(processes=settings["CPUs"])
            results = [pool.apply_async(self.cutadapt_PE, \
                              args = (settings,library)) \
                              for library in sorted(settings["libraries"])]
            pool.close()
            pool.join()
            for r in results:
                r.get()
        else:
            pool = mp.Pool(processes=settings["CPUs"])
            results = [pool.apply_async(self.cutadapt_SE, \
                              args = (settings,library)) \
                              for library in sorted(settings["libraries"])]
            pool.close()
            pool.join()
            for r in results:
                r.get()
        self.test_errors(settings)

    def make_folder(self,settings):
        '''
        make output folder for task if needed
        '''
        if not os.path.exists(os.path.join(settings["--output"],"trim")):
            os.makedirs(os.path.join(settings["--output"],"trim"))

        if not os.path.exists(os.path.join(settings["--output"],"trim","triminfo")):
            os.makedirs(os.path.join(settings["--output"],"trim","triminfo"))
        if not os.path.exists(os.path.join(settings["--output"],"trim","discard")):
            os.makedirs(os.path.join(settings["--output"],"trim","discard"))
            
    def cutadapt_PE(self,settings,library):
        '''
        Run cutadapt with paired end settings
        '''
        print("\t"+library)
        input_file1 = os.path.join(settings["--input"],library+\
                                  settings["trim"]["trim_input_file_suffix_for"])
        input_file2 = os.path.join(settings["--input"],library+\
                                  settings["trim"]["trim_input_file_suffix_rev"])
        output_file = os.path.join(settings["--output"],"trim",library)
        short_folder = os.path.join(settings["--output"],"trim",\
                                    "discard",library)
        trim_info_file = os.path.join(settings["--output"],"trim",\
                                      "triminfo",library+"_triminfo.log")
        trim_error_file = os.path.join(settings["--output"],"trim",\
                                      "triminfo",library+"_trimerror.log")
        cutadapt_command = (
                    settings["trim"]["trim_call"],
                    "-a", "file:"+settings["trim"]["trim_adapter_for"],
                    "-A", "file:"+settings["trim"]["trim_adapter_rev"],
                    "--too-short-output", short_folder+"_1_short.fq",
                    "--too-short-paired-output", short_folder+"_2_short.fq",
                    "--minimum-length", str(settings["min_length"]),
                    "-q", str(settings["trim"]["trim_min_qual"]),
                    "--quality-base", str(settings["trim"]["trim_quality_base"]),
                    "-O", str(settings["trim"]["trim_overlap"]),
                    "-j", str(settings["trim"]["trim_threads"]),
                    "-o", output_file+"_1_trim.fq",
                    "-p", output_file+"_2_trim.fq",
                    input_file1,input_file2,
                    "2>&1", ">", trim_info_file, "|", "tee", trim_error_file
                    )
        os.system(" ".join(cutadapt_command))
    
    def cutadapt_SE(self,settings,library):
        '''
        Run cutadapt with paired end settings
        '''
        print("\t"+library)
        input_file = os.path.join(settings["--input"],library+\
                                  settings["trim"]["trim_input_file_suffix_SE"])
        output_file = os.path.join(settings["--output"],"trim",library)
        short_folder = os.path.join(settings["--output"],"trim",\
                                    "discard",library)
        trim_info_file = os.path.join(settings["--output"],"trim",\
                                      "triminfo",library+"_triminfo.txt")
        trim_error_file = os.path.join(settings["--output"],"trim",\
                                      "triminfo",library+"_trimerror.log")
        cutadapt_command = (
                    settings["trim"]["trim_call"],
                    "-a", "file:"+settings["trim"]["trim_adapter_for"],
                    "--too-short-output", short_folder+"_short.fq",
                    "--untrimmed-output", short_folder+"_untrimmed.fq",
                    "--minimum-length", str(settings["min_length"]),
                    "-q", str(settings["trim"]["trim_min_qual"]),
                    "--quality-base", str(settings["trim"]["trim_quality_base"]),
                    "-j", str(settings["trim"]["trim_threads"]),
                    "-O", str(max(settings["trim"]["trim_overlap"],3)),
                    "-o", output_file+"_trim.fq", input_file,
                    "2>&1", ">", trim_info_file, "|", "tee", trim_error_file
                    )
        os.system(" ".join(cutadapt_command))

    def test_errors(self,settings):
        '''
        Tests are any error files containing a error raport
        '''
        info_folder = os.path.join(settings["--output"],"trim",\
                                          "triminfo")
        if len(os.listdir(info_folder)) == 0:
            sys.exit('Task "trim" incomplete! Infofiles missing in folder ' + info_folder)
        for file in os.listdir(info_folder):
            if file.endswith("_trimerror.log"):
                if os.stat(os.path.join(info_folder,file)).st_size != 0:
                    sys.exit('Error in task "trim", check infofile: '+file)
        
        


        
