# coding=utf-8

"""
"""
import os
import unittest
import shutil
from common_wrangler.add2dict import main
from common_wrangler.common import (capture_stdout, capture_stderr, silent_remove, diff_lines)
import logging

# logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
DISABLE_REMOVE = logger.isEnabledFor(logging.DEBUG)

__author__ = 'hmayes'


# Directories #

DATA_DIR = os.path.join(os.path.dirname(__file__), 'test_data')
SUB_DATA_DIR = os.path.join(DATA_DIR, 'add2dict')

# Files #

CUST_DICT = os.path.join(SUB_DATA_DIR, 'cust.dic')
TEMP_CUST_DICT = os.path.join(SUB_DATA_DIR, 'temp_cust.dic')
GOOD_CUST_DICT1 = os.path.join(SUB_DATA_DIR, 'cust1_good.dic')


class TestAdd2DictNoOutput(unittest.TestCase):
    def testHelp(self):
        test_input = ['-h']
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertFalse(output)
        with capture_stdout(main, test_input) as output:
            self.assertTrue("optional arguments" in output)

    def testIOErrorDict(self):
        test_input = ["oligomer", "-d", "ghost.dic"]
        # main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("No such file or directory:" in output)

    def testMissingArg(self):
        test_input = ["-d", TEMP_CUST_DICT]
        # main(test_input)
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("following arguments are required" in output)


class TestAdd2Dict(unittest.TestCase):
    def testNormalUse(self):
        shutil.copyfile(CUST_DICT, TEMP_CUST_DICT)
        test_input = ["oligomer", "-d", TEMP_CUST_DICT]
        try:
            main(test_input)
        finally:
            silent_remove(TEMP_CUST_DICT, DISABLE_REMOVE)
            pass
