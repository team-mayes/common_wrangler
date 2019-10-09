# coding=utf-8

"""
"""
import os
import unittest

import shutil

from common_wrangler.rename_files import main
from common_wrangler.common import (capture_stdout, capture_stderr, silent_remove)
import logging

# logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
DISABLE_REMOVE = logger.isEnabledFor(logging.DEBUG)

__author__ = 'hmayes'


# Directories #

DATA_DIR = os.path.join(os.path.dirname(__file__), 'test_data')
SUB_DATA_DIR = os.path.join(DATA_DIR, 'rename_files')
SUB_SUB_DATA_DIR = os.path.join(SUB_DATA_DIR, 'sub_dir')

# Files #

SMALL_FILE = os.path.join(SUB_DATA_DIR, 'small_file.txt')


# test data #
TEST_FILE_NAMES = ['has space.txt', 'has two spaces.txt', 'now!exclaim.txt']
REPLACED_FILE_NAMES1 = ['hasspace.txt', 'hastwospaces.txt', 'now!exclaim.txt']
REPLACED_FILE_NAMES2 = ['has space.txt', 'has two spaces.txt', 'now_exclaim.txt']

# REPLACED_FILE_NAMES3 = ['has_space.txt', 'has_two_spaces.txt', 'now!exclaim.txt']


def make_files(fname_list, file_dir):
    """
    Create files fresh, because will be moved when program runs
    @param fname_list: list of file names without directory name
    @param file_dir: name of the directory where the files should be created
    """
    for fname in fname_list:
        new_file = os.path.join(file_dir, fname)
        shutil.copyfile(SMALL_FILE, new_file)


def get_abs_path_name(fname_list, abs_dir):
    """
    Create files fresh, because will be moved when program runs
    @param fname_list: list of file names without directory name
    @param abs_dir: absolute directory name
    @return full_name_list: a list of file names with the specified absolute directory
    """
    full_name_list = []
    for fname in fname_list:
        full_name_list.append(os.path.join(abs_dir, fname))
    return full_name_list


def count_files(fname_list):
    """
    Counts how many files in list exist
    @param fname_list: list of file names
    @return num_existing_files: a list of file names with the specified absolute directory
    """
    num_existing_files = 0
    for fname in fname_list:
        if os.path.isfile(fname):
            num_existing_files += 1
    return num_existing_files


class TestRenameNoOutput(unittest.TestCase):
    def testHelp(self):
        test_input = ['-h']
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertFalse(output)
        with capture_stdout(main, test_input) as output:
            self.assertTrue("optional arguments" in output)

    def testInvalidArg(self):
        test_input = ['-@']
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("unrecognized arguments" in output)

    def testNoFilesRenamed(self):
        test_input = []
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stdout(main, test_input) as output:
            self.assertTrue("Found and renamed 0 files" in output)


class TestRename(unittest.TestCase):
    def testDefaultPatterns(self):
        # test with a sub dir
        initial_fnames = get_abs_path_name(TEST_FILE_NAMES, SUB_DATA_DIR) + get_abs_path_name(TEST_FILE_NAMES,
                                                                                              SUB_SUB_DATA_DIR)
        expected_fnames = get_abs_path_name(REPLACED_FILE_NAMES1,
                                            SUB_DATA_DIR) + get_abs_path_name(REPLACED_FILE_NAMES1, SUB_SUB_DATA_DIR)
        for fname in expected_fnames + initial_fnames:
            silent_remove(fname, disable=DISABLE_REMOVE)
        make_files(TEST_FILE_NAMES, SUB_DATA_DIR)
        make_files(TEST_FILE_NAMES, SUB_SUB_DATA_DIR)
        test_input = ["-d", SUB_DATA_DIR]
        try:
            if logger.isEnabledFor(logging.DEBUG):
                main(test_input)
                # need to make again after test run for capturing std out
                make_files(TEST_FILE_NAMES, SUB_DATA_DIR)
                make_files(TEST_FILE_NAMES, SUB_SUB_DATA_DIR)
            # main(test_input)
            with capture_stdout(main, test_input) as output:
                self.assertTrue("Found and renamed 4 files" in output)
            self.assertTrue(count_files(initial_fnames) == 2)
            self.assertTrue(count_files(expected_fnames) == 6)
        finally:
            for fname in expected_fnames:
                silent_remove(fname, disable=DISABLE_REMOVE)
            pass

    def testDefaultPatternsOnlyOneDeep(self):
        # test with a sub dir
        initial_fnames = get_abs_path_name(TEST_FILE_NAMES, SUB_DATA_DIR) + get_abs_path_name(TEST_FILE_NAMES,
                                                                                              SUB_SUB_DATA_DIR)
        expected_fnames = get_abs_path_name(REPLACED_FILE_NAMES1,
                                            SUB_DATA_DIR) + get_abs_path_name(TEST_FILE_NAMES, SUB_SUB_DATA_DIR)
        for fname in expected_fnames + initial_fnames:
            silent_remove(fname, disable=DISABLE_REMOVE)
        make_files(TEST_FILE_NAMES, SUB_DATA_DIR)
        make_files(TEST_FILE_NAMES, SUB_SUB_DATA_DIR)
        test_input = ["-d", SUB_DATA_DIR, "-o"]
        try:
            if logger.isEnabledFor(logging.DEBUG):
                main(test_input)
                # need to make again after test run for capturing std out
                make_files(TEST_FILE_NAMES, SUB_DATA_DIR)
                make_files(TEST_FILE_NAMES, SUB_SUB_DATA_DIR)
            with capture_stdout(main, test_input) as output:
                self.assertTrue("Found and renamed 2 files" in output)
            self.assertTrue(count_files(initial_fnames) == 4)
            self.assertTrue(count_files(expected_fnames) == 6)
        finally:
            for fname in expected_fnames:
                silent_remove(fname, disable=DISABLE_REMOVE)
            pass

    def testAltPattern(self):
        initial_fnames = get_abs_path_name(TEST_FILE_NAMES, SUB_DATA_DIR)
        expected_fnames = get_abs_path_name(REPLACED_FILE_NAMES2, SUB_DATA_DIR)
        for fname in expected_fnames + initial_fnames:
            silent_remove(fname, disable=DISABLE_REMOVE)
        make_files(TEST_FILE_NAMES, SUB_DATA_DIR)
        test_input = ["-d", SUB_DATA_DIR, "-p", "!", "-n", "_"]
        try:
            if logger.isEnabledFor(logging.DEBUG):
                main(test_input)
                # need to make again for capturing std out
                make_files(TEST_FILE_NAMES, SUB_DATA_DIR)
            with capture_stdout(main, test_input) as output:
                self.assertTrue("Found and renamed 1 files" in output)
            self.assertTrue(count_files(initial_fnames), 1)
            self.assertTrue(count_files(expected_fnames), 3)
        finally:
            for fname in expected_fnames:
                silent_remove(fname, disable=DISABLE_REMOVE)
