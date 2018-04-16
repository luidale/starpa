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


    def test_01(self):
        input_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","fq")))
        output_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output")))
        args = docopt(doc, ["-s","trim","-e","trim","-c",config_file,\
                            "-i",input_folder,"-o","tests/data/output"])
        #change trim call in Travis
        if os.environ.get('TRAVIS') == 'true':
            args["cutaadapt"]["trim_call"] = "cutadapt"
        starpa.main(args)
        shutil.rmtree("tests/data/output")
        #self.assertTrue(run)

    def test_02(self):
        input_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","fq")))
        output_folder = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           os.path.join("data","output")))
        args = docopt(doc, ["-s","trim","-e","quantify","-c",config_file,\
                            "-i",input_folder,"-o","tests/data/output"])
        starpa.main(args)
        shutil.rmtree("tests/data/output")
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
