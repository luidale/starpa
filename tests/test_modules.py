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


class TestStarpa(unittest.TestCase):

    def setUp(self):
        pass


##    def test_01(self):
##        input_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
##                           os.path.join("data","fq")))
##        output_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
##                           os.path.join("data","output")))
##        args = docopt(doc, ["-s","trim","-e","trim","-c",config_file,\
##                            "-i",input_folder,"-o","tests/data/output"])
##        print(args)
##        #change trim call in Travis
####        if os.environ.get('TRAVIS') == 'true':
####            args["cutadapt"]["trim_call"] = "cutadapt2"
##        starpa.main(args)
##        shutil.rmtree("tests/data/output")
##        #self.assertTrue(run)

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
        if os.path.isdir(output_folder):
            print(output_folder)
            if os.path.isdir(os.path.join(output_folder,"trim")):
                print(os.path.join(output_folder,"trim"))
            else:
                print("A")
        else:
            print("B")
            

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
        input_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","fq")))
        output_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output")))
        args = docopt(doc, ["-s","trim","-e","quantify", "-c", config_file,\
                            "-i",input_folder,"-o",output_folder])
        starpa.main(args)
        shutil.rmtree(output_folder)
        #self.assertTrue(run)

    def test_01a(self):
        args = docopt(doc, ["--version"])
        #print(args[1])
        #run = starpa.main(["--start","align"])
        self.assertEqual(args["--version"],True)


def main():
    unittest.main()

if __name__ == '__main__':
    main()
