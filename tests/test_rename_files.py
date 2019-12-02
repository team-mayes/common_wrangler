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
SUB_SUB_DIR = os.path.join(SUB_DATA_DIR, 'sub_dir')

# Files #

SMALL_FILE = os.path.join(SUB_DATA_DIR, 'small_file.txt')


# test data #
TEST_FNAMES = ['has space.txt', 'has two spaces.txt', 'now!exclaim.txt', 'WOW.csv']
MORE_FNAMES = ['wow_wow.CSV', 'WOW-WOW-WOW.CSV']
REPLACED_FNAMES1 = ['hasspace.txt', 'hastwospaces.txt', 'now!exclaim.txt', 'WOW.csv']
REPLACED_FNAMES2 = ['has space.txt', 'has two spaces.txt', 'now_exclaim.txt', 'WOW.csv']
REPLACED_FNAMES3 = ['has space.txt', 'has two spaces.txt', 'now!exclaim.txt', 'wow.csv', 'wow_wow.csv',
                    'wow-wow-wow.csv']
CLEAN_UP_FNAMES = set(TEST_FNAMES + REPLACED_FNAMES1 + REPLACED_FNAMES2 + REPLACED_FNAMES3)
# REPLACED_FILE_NAMES3 = ['has_space.txt', 'has_two_spaces.txt', 'now!exclaim.txt']


def get_rename_files_test_list_to_clean():
    to_clean_fnames = get_abs_path(CLEAN_UP_FNAMES, SUB_DATA_DIR) + get_abs_path(CLEAN_UP_FNAMES, SUB_SUB_DIR)
    return to_clean_fnames + [SUB_SUB_DIR]


def make_files(fname_list, file_dir):
    """
    Create files fresh, because will be moved when program runs
    :param fname_list: list of file names without directory name
    :param file_dir: name of the directory where the files should be created
    :return: list of locations of created files
    """
    initial_fnames = []
    for fname in fname_list:
        new_file_name = os.path.join(file_dir, fname)
        shutil.copyfile(SMALL_FILE, new_file_name)
        initial_fnames.append(new_file_name)
    return initial_fnames


def get_abs_path(fname_list, abs_dir):
    """
    Create files fresh, because will be moved when program runs
    :param fname_list: list of file names without directory name
    :param abs_dir: absolute directory name
    :return: full_name_list: a list of file names with the specified absolute directory
    """
    full_name_list = []
    for fname in fname_list:
        base_name = os.path.basename(fname)
        full_name_list.append(os.path.join(abs_dir, base_name))
    return full_name_list


def count_files(fname_list):
    """
    Counts how many files in list exist
    :param fname_list: list of file names
    :return: num_existing_files: a list of file names with the specified absolute directory
    """
    num_existing_files = 0
    for fname in fname_list:
        if os.path.isfile(fname):
            num_existing_files += 1
    return num_existing_files


def clean_then_make_files(list_to_clean, fname_to_make_list, dir_list):
    """
    Put this repeated task in a method
    :param list_to_clean: list with file and/or directory names
    :param fname_to_make_list: list of file names
    :param dir_list: list of one or more names of directories where files should be created
    :return: initial_fnames: list of locations of created files
    """
    # start clean
    for f_or_dir_name in list_to_clean:
        silent_remove(f_or_dir_name)
    # make files
    initial_fnames = []
    for dir_name in dir_list:
        if dir_name:
            if not os.path.exists(dir_name):
                os.makedirs(dir_name)
            initial_fnames += make_files(fname_to_make_list, dir_name)
    return initial_fnames


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
        # since this would fail if didn't start from clean state, make sure no remain files from any other test
        list_to_clean = get_rename_files_test_list_to_clean()
        for f_or_dir_name in list_to_clean:
            silent_remove(f_or_dir_name)

        test_input = []
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stdout(main, test_input) as output:
            self.assertTrue("Found and renamed 0 files" in output)


class TestRename(unittest.TestCase):
    def testDefaultPatterns(self):
        # test with a sub dir
        list_to_clean = get_rename_files_test_list_to_clean()
        initial_fnames = clean_then_make_files(list_to_clean, TEST_FNAMES, [SUB_DATA_DIR, SUB_SUB_DIR])
        expected_fnames = get_abs_path(REPLACED_FNAMES1, SUB_DATA_DIR) + get_abs_path(REPLACED_FNAMES1, SUB_SUB_DIR)
        test_input = ["-d", SUB_DATA_DIR]
        try:
            if logger.isEnabledFor(logging.DEBUG):
                main(test_input)
                # reset for next run, first cleaning and the recreating files
                clean_then_make_files(list_to_clean, TEST_FNAMES, [SUB_DATA_DIR, SUB_SUB_DIR])
            # main(test_input)
            with capture_stdout(main, test_input) as output:
                self.assertTrue("Found and renamed 4 files" in output)
            self.assertTrue(count_files(initial_fnames) == 4)
            self.assertTrue(count_files(expected_fnames) == 8)
        finally:
            for f_or_dir_name in list_to_clean:
                silent_remove(f_or_dir_name, DISABLE_REMOVE)
            pass

    def testDefaultPatternsOnlyOneDeep(self):
        # test with a sub dir
        list_to_clean = get_rename_files_test_list_to_clean()
        initial_fnames = clean_then_make_files(list_to_clean, TEST_FNAMES, [SUB_DATA_DIR, SUB_SUB_DIR])
        expected_fnames = get_abs_path(REPLACED_FNAMES1, SUB_DATA_DIR) + get_abs_path(TEST_FNAMES, SUB_SUB_DIR)
        test_input = ["-d", SUB_DATA_DIR, "-o"]
        try:
            if logger.isEnabledFor(logging.DEBUG):
                main(test_input)
                # need to make again after test run for capturing std out
                clean_then_make_files(list_to_clean, TEST_FNAMES, [SUB_DATA_DIR, SUB_SUB_DIR])
            with capture_stdout(main, test_input) as output:
                self.assertTrue("Found and renamed 2 files" in output)
            self.assertTrue(count_files(initial_fnames) == 6)
            self.assertTrue(count_files(expected_fnames) == 8)
        finally:
            for f_or_dir_name in list_to_clean:
                silent_remove(f_or_dir_name, DISABLE_REMOVE)
            pass

    def testAltPattern(self):
        list_to_clean = get_rename_files_test_list_to_clean()
        initial_fnames = clean_then_make_files(list_to_clean, TEST_FNAMES, [SUB_DATA_DIR])
        expected_fnames = get_abs_path(REPLACED_FNAMES2, SUB_DATA_DIR)
        test_input = ["-d", SUB_DATA_DIR, "-p", "!", "-n", "_"]
        try:
            if logger.isEnabledFor(logging.DEBUG):
                main(test_input)
                # need to make again for capturing std out
                clean_then_make_files(list_to_clean, TEST_FNAMES, [SUB_DATA_DIR])
            with capture_stdout(main, test_input) as output:
                self.assertTrue("Found and renamed 1 files" in output)
            self.assertTrue(count_files(initial_fnames), 3)
            self.assertTrue(count_files(expected_fnames), 4)
        finally:
            for f_or_dir_name in list_to_clean:
                silent_remove(f_or_dir_name, DISABLE_REMOVE)

    def testMakeLowerCase(self):
        list_to_clean = get_rename_files_test_list_to_clean()
        initial_fnames = clean_then_make_files(list_to_clean, TEST_FNAMES + MORE_FNAMES, [SUB_DATA_DIR])
        expected_fnames = get_abs_path(REPLACED_FNAMES3, SUB_DATA_DIR)
        test_input = ["-d", SUB_DATA_DIR, "-l"]
        try:
            if logger.isEnabledFor(logging.DEBUG):
                main(test_input)
                # need to make again for capturing std out
                clean_then_make_files(list_to_clean, TEST_FNAMES, [SUB_DATA_DIR])
            with capture_stdout(main, test_input) as output:
                self.assertTrue("Found and renamed 3 files" in output)
            self.assertTrue(count_files(initial_fnames), 6)
            self.assertTrue(count_files(expected_fnames), 6)
        finally:
            for f_or_dir_name in list_to_clean:
                silent_remove(f_or_dir_name, DISABLE_REMOVE)

    def testMakeLowercaseTestingMode(self):
        list_to_clean = get_rename_files_test_list_to_clean()
        initial_fnames = clean_then_make_files(list_to_clean, TEST_FNAMES + MORE_FNAMES, [SUB_DATA_DIR])
        expected_fnames = initial_fnames
        test_input = ["-d", SUB_DATA_DIR, "-l", "-t"]
        try:
            if logger.isEnabledFor(logging.DEBUG):
                main(test_input)
                # need to make again for capturing std out
                clean_then_make_files(list_to_clean, TEST_FNAMES, [SUB_DATA_DIR])
            with capture_stdout(main, test_input) as output:
                self.assertTrue("would rename 3 files" in output)
            self.assertTrue(count_files(initial_fnames), 6)
            self.assertTrue(count_files(expected_fnames), 6)
        finally:
            for f_or_dir_name in list_to_clean:
                silent_remove(f_or_dir_name, DISABLE_REMOVE)
