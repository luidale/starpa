#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Tests should be run from main folder
'''

import sys
import os
import shutil
#allows to import from the src folder
sys.path[0] = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("..","src")))
import starpa.__main__ as starpa
import unittest
from docopt import docopt

doc = starpa.__doc__
config_file = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","config.txt")))
config_file_SE = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","config_SE.txt")))
config_file_sens = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","config_sens.txt")))
config_file_sens_SE = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","config_sens_SE.txt")))

class TestStarpa(unittest.TestCase):

    def setUp(self):
        pass

    def test_01_trim(self):
        if '__pypy__' in sys.builtin_module_names:
            return
        input_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","fq")))
        output_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output")))
        args = docopt(doc, ["-s","trim","-e","trim","-c",config_file,\
                            "-i",input_folder,"-o",output_folder])
        starpa.main(args)

    def test_02_align(self):

        output_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output")))
        if '__pypy__' in sys.builtin_module_names:
            #in pypy the cutadapt does not work, prepared input is used
            input_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","trim_output")))

        else:
            input_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output","trim")))
            
        args = docopt(doc, ["-s","align","-e","align","-c",config_file,\
                            "-i",input_folder,"-o",output_folder])
        starpa.main(args)

    def test_02sens_align(self):

        output_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output")))
        if '__pypy__' in sys.builtin_module_names:
            #in pypy the cutadapt does not work, prepared input is used
            input_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","trim_output")))

        else:
            input_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output","trim")))
            
        args = docopt(doc, ["-s","align","-e","align","-c",config_file_sens,\
                            "-i",input_folder,"-o",output_folder])
        starpa.main(args)


    def test_03_sort(self):
        input_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output","align")))
        output_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output")))
        args = docopt(doc, ["-s","sam_sort","-e","sam_sort","-c",config_file,\
                            "-i",input_folder,"-o",output_folder])
        starpa.main(args)

    def test_04_pseudoSE(self):
        input_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output","sam_sort")))
        output_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output")))
        args = docopt(doc, ["-s","pseudoSE","-e","pseudoSE","-c",config_file,\
                            "-i",input_folder,"-o",output_folder])
        starpa.main(args)

    def test_05_identify(self):
        input_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output","pseudoSE")))
        output_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output")))
        args = docopt(doc, ["-s","identify","-e","identify","-c",config_file,\
                            "-i",input_folder,"-o",output_folder])
        starpa.main(args)

    def test_06_cluster(self):
        input_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output","identify")))
        output_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output")))
        args = docopt(doc, ["-s","cluster","-e","cluster","-c",config_file,\
                            "-i",input_folder,"-o",output_folder])
        starpa.main(args)

    def test_07_quantify(self):
        input_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output","cluster")))
        output_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output")))
        args = docopt(doc, ["-s","quantify","-e","quantify","-c",config_file,\
                            "-i",input_folder,"-o",output_folder])
        starpa.main(args)
        shutil.rmtree(output_folder)
        
    def test_08_full(self):
        if '__pypy__' in sys.builtin_module_names:
            #in pypy the cutadapt does not work, prepared input is used
            input_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","trim_output")))
            start_task = "align"
        else:
            input_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","fq")))
            start_task = "trim"
        output_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output")))
        args = docopt(doc, ["-s",start_task,"-e","quantify", "-c", config_file,\
                            "-i",input_folder,"-o",output_folder])
        starpa.main(args)
        shutil.rmtree(output_folder)
        #self.assertTrue(run)

    def test_09_SE_trim(self):
        if '__pypy__' in sys.builtin_module_names:
            return
        input_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","fq")))
        output_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output")))
        args = docopt(doc, ["-s","trim","-e","trim","-c",config_file_SE,\
                            "-i",input_folder,"-o",output_folder])
        starpa.main(args)

    def test_10_SE_align(self):

        output_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output")))
        if '__pypy__' in sys.builtin_module_names:
            #in pypy the cutadapt does not work, prepared input is used
            input_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","trim_SE_output")))

        else:
            input_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output","trim")))
            
        args = docopt(doc, ["-s","align","-e","align","-c",config_file_SE,\
                            "-i",input_folder,"-o",output_folder])
        starpa.main(args)

    def test_10_SE_align(self):

        output_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output")))
        if '__pypy__' in sys.builtin_module_names:
            #in pypy the cutadapt does not work, prepared input is used
            input_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","trim_SE_output")))

        else:
            input_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output","trim")))
            
        args = docopt(doc, ["-s","align","-e","align","-c",config_file_sens_SE,\
                            "-i",input_folder,"-o",output_folder])
        starpa.main(args)
        
    def test_11_SE_sort(self):
        input_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output","align")))
        output_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output")))
        args = docopt(doc, ["-s","sam_sort","-e","sam_sort","-c",config_file_SE,\
                            "-i",input_folder,"-o",output_folder])
        starpa.main(args)

    def test_12_SE_pseudoSE(self):
        input_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output","sam_sort")))
        output_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output")))
        args = docopt(doc, ["-s","pseudoSE","-e","pseudoSE","-c",config_file_SE,\
                            "-i",input_folder,"-o",output_folder])
        starpa.main(args)

    def test_01a(self):
        args = docopt(doc, ["--version"])
        #print(args[1])
        #run = starpa.main(["--start","align"])
        self.assertEqual(args["--version"],True)


def main():
    unittest.main()

if __name__ == '__main__':
    main()
