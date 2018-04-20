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

    def test_trim(self):
        input_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","fq")))
        output_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output")))
        args = docopt(doc, ["-s","trim","-e","trim","-c",config_file,\
                            "-i",input_folder,"-o",output_folder])
        starpa.main(args)

    def test_align(self):
        input_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output","trim")))
        output_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output")))
        args = docopt(doc, ["-s","align","-e","align","-c",config_file,\
                            "-i",input_folder,"-o",output_folder])
        starpa.main(args)

    def test_sort(self):
        input_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output","align")))
        output_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output")))
        args = docopt(doc, ["-s","sam_sort","-e","sam_sort","-c",config_file,\
                            "-i",input_folder,"-o",output_folder])
        starpa.main(args)

    def test_pseudoSE(self):
        input_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output","sam_sort")))
        output_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output")))
        args = docopt(doc, ["-s","pseudoSE","-e","pseudoSE","-c",config_file,\
                            "-i",input_folder,"-o",output_folder])
        starpa.main(args)

##    def test_identify(self):
##        input_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
##                           os.path.join("data","output","pseudoSE")))
##        output_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
##                           os.path.join("data","output")))
##        args = docopt(doc, ["-s","identify","-e","identify","-c",config_file,\
##                            "-i",input_folder,"-o","tests/data/output"])
##        starpa.main(args)
##
##    def test_cluster(self):
##        input_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
##                           os.path.join("data","output","identify")))
##        output_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
##                           os.path.join("data","output")))
##        args = docopt(doc, ["-s","cluster","-e","cluster","-c",config_file,\
##                            "-i",input_folder,"-o","tests/data/output"])
##        starpa.main(args)
##
##    def test_quantify(self):
##        input_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
##                           os.path.join("data","output","cluster")))
##        output_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
##                           os.path.join("data","output")))
##        args = docopt(doc, ["-s","quantify","-e","quantify","-c",config_file,\
##                            "-i",input_folder,"-o","tests/data/output"])
##        starpa.main(args)
##        shutil.rmtree("tests/data/output")
##        
##    def test_full(self):
##        input_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
##                           os.path.join("data","fq")))
##        output_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
##                           os.path.join("data","output")))
##        args = docopt(doc, ["-s","trim","-e","quantify", "-c", config_file,\
##                            "-i",input_folder,"-o","tests/data/output"])
##        starpa.main(args)
##        shutil.rmtree("tests/data/output")
##        #self.assertTrue(run)

    def test_01a(self):
        args = docopt(doc, ["--version"])
        #print(args[1])
        #run = starpa.main(["--start","align"])
        self.assertEqual(args["--version"],True)


def main():
    unittest.main()

if __name__ == '__main__':
    main()
