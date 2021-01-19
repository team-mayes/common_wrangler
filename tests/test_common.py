#!/usr/bin/env python3
"""
Unit and regression test for the common_wrangler package.
"""

# Import package, test suite, and other packages as needed
import numpy as np
import os
import shutil
import tempfile
import unittest
from argparse import Namespace
from common_wrangler.common import (NUM_ATOMS, MAIN_SEC, SEC_ATOMS, SEC_HEAD, SEC_TAIL, ATOM_COORDS, ATOM_TYPE,
                                    COLOR_SEQUENCE, InvalidDataError, NotFoundError, TemplateNotReadableError,
                                    find_files_by_dir, read_csv, get_fname_root, write_csv, str_to_bool,
                                    read_csv_header, fmt_row_data, calc_k, diff_lines, create_out_fname, dequote, quote,
                                    conv_raw_val, pbc_calc_vector, pbc_vector_avg, unit_vector, vec_angle, vec_dihedral,
                                    check_for_files, make_dir, silent_remove, list_to_file, list_to_csv,
                                    longest_common_substring, capture_stdout, print_csv_stdout, read_tpl,
                                    file_rows_to_list, round_to_12th_decimal, single_quote, calc_dist,
                                    np_float_array_from_file, capture_stderr, round_sig_figs, process_cfg, parse_stoich,
                                    natural_keys, str_to_file, read_csv_to_list, read_csv_dict, round_to_fraction,
                                    make_fig, read_json, process_pdb_file, assign_color, conv_str_to_func, unique_list,
                                    overwrite_config_vals)
import logging

try:
    # noinspection PyCompatibility
    from ConfigParser import ConfigParser, MissingSectionHeaderError
except ImportError:
    # noinspection PyCompatibility
    from configparser import ConfigParser, MissingSectionHeaderError


__author__ = 'mayes'

# logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
DISABLE_REMOVE = logger.isEnabledFor(logging.DEBUG)

# Constants #
DATA_DIR = os.path.join(os.path.dirname(__file__), 'test_data')
SUB_DATA_DIR = os.path.join(DATA_DIR, 'common')
FES_DIR = os.path.join(SUB_DATA_DIR, 'fes_out')
NEW_DIR = os.path.join(SUB_DATA_DIR, 'new_dir')

ELEM_DICT_FILE = os.path.join(SUB_DATA_DIR, 'element_dict.csv')
ATOM_DICT_FILE = os.path.join(SUB_DATA_DIR, 'atom_reorder.csv')
GOOD_ATOM_DICT = {1: 20, 2: 21, 3: 22, 4: 23, 5: 24, 6: 25, 7: 26, 8: 27, 9: 2, 10: 1, 11: 3, 12: 4, 13: 5, 14: 6,
                  15: 7, 16: 8, 17: 9, 18: 10, 19: 11, 20: 12, 21: 13, 22: 14, 23: 15, 24: 16, 25: 17, 26: 18, 27: 19}

CSV_FILE = os.path.join(DATA_DIR, SUB_DATA_DIR, 'rad_PMF_last2ns3_1.txt')
ALT_CSV_FILE = os.path.join(DATA_DIR, SUB_DATA_DIR, 'rad_PMF_last2ns3_1_alt.txt')
FRENG_TYPES = [float, str]

ONE_KEY_INI = os.path.join(SUB_DATA_DIR, 'one_key_config.ini')

ORIG_WHAM_ROOT = "PMF_last2ns3_1"
ORIG_WHAM_FNAME = ORIG_WHAM_ROOT + ".txt"
ORIG_WHAM_PATH = os.path.join(DATA_DIR, ORIG_WHAM_FNAME)
SHORT_WHAM_PATH = os.path.join(DATA_DIR, ORIG_WHAM_FNAME)
EMPTY_CSV = os.path.join(SUB_DATA_DIR, 'empty.csv')

FILE_LIST = os.path.join(SUB_DATA_DIR, 'file_list.txt')
FILE_LIST_W_MISSING_FILE = os.path.join(SUB_DATA_DIR, 'file_list_with_ghost.txt')

TEST_DIR_SEARCH_FILE = os.path.join(NEW_DIR, 'business_file.csv')
TEST_DIR_SEARCH_FILE2 = os.path.join(NEW_DIR, 'business_file.unique')

BOX_SIZES_NO_HEADER_SPACE_SEP = os.path.join(SUB_DATA_DIR, 'qm_box_sizes.txt')
BOX_SIZES_HEADER_SPACE_SEP = os.path.join(SUB_DATA_DIR, 'qm_box_sizes_header.txt')
BOX_SIZES_NO_HEADER_COMMA_SEP = os.path.join(SUB_DATA_DIR, 'qm_box_sizes.csv')
BOX_SIZES_HEADER_COMMA_SEP = os.path.join(SUB_DATA_DIR, 'qm_box_sizes_header.csv')

GOOD_BOX_NDARRAY = np.asarray([[11.89100027, 15.36799955, 18.10500002], [11.375, 14.99599981, 17.98800039],
                               [11.10300016, 15.60500002, 18.31400013], [10., 14.995, 10.98800039]])
GOOD_BOX_NDARRAY_ROW_NAN = np.asarray([[np.nan, np.nan, np.nan],
                                       [11.89100027, 15.36799955, 18.10500002], [11.375, 14.99599981, 17.98800039],
                                       [11.10300016, 15.60500002, 18.31400013], [10., 14.995, 10.98800039]])

FLOAT_AND_NON = os.path.join(SUB_DATA_DIR, 'msm_sum_output.csv')

TEST_PLOT_FNAME = os.path.join(SUB_DATA_DIR, 'sample_plot.png')

OUT_PFX = 'rad_'

# Data #

CSV_HEADER = ['coord', 'free_energy', 'corr']
GHOST = 'ghost'
NUM_COLORS = 11

# Output files #

DIFF_LINES_BASE_FILE = os.path.join(SUB_DATA_DIR, 'diff_lines_base_file.csv')
DIFF_LINES_PREC_DIFF = os.path.join(SUB_DATA_DIR, 'diff_lines_prec_diff.csv')
DIFF_LINES_ONE_VAL_DIFF = os.path.join(SUB_DATA_DIR, 'diff_lines_one_val_diff.csv')
DIFF_LINES_MISS_VAL = os.path.join(SUB_DATA_DIR, 'diff_lines_miss_val.csv')
MISS_LINES_MISS_LINE = os.path.join(SUB_DATA_DIR, 'diff_lines_miss_line.csv')
DIFF_LINES_ONE_NAN = os.path.join(SUB_DATA_DIR, 'diff_lines_one_nan.csv')
DIFF_LINES_ONE_NAN_PREC_DIFF = os.path.join(SUB_DATA_DIR, 'diff_lines_one_nan.csv')
DIFF_LINES_STR_DIFF = os.path.join(SUB_DATA_DIR, 'diff_lines_str_diff.csv')
DIFF_LINES_DIFF_ORDER = os.path.join(SUB_DATA_DIR, 'diff_lines_diff_order.csv')

DIFF_LINES_SCI_FILE = os.path.join(SUB_DATA_DIR, 'cv_analysis_quat.log')
DIFF_LINES_ALT_SCI_FILE = os.path.join(SUB_DATA_DIR, 'cv_analysis_quat_good.log')

LIST_OUT = os.path.join(SUB_DATA_DIR, "temp.txt")
GOOD_LIST_OUT = os.path.join(SUB_DATA_DIR, "good_list.txt")
GOOD_FORMAT_LIST_OUT = os.path.join(SUB_DATA_DIR, "good_format_list.txt")

DEF_FILE_PAT = 'fes*.out'

IMPROP_SEC = os.path.join(SUB_DATA_DIR, 'glue_improp.data')
IMPROP_SEC_ALT = os.path.join(SUB_DATA_DIR, 'glue_improp_diff_ord.data')

# To test PBC math
PBC_BOX = np.full(3, 24.25)
A_VEC = [3.732, -1.803, -1.523]
B_VEC = [4.117, 0.135, -2.518]
GOOD_A_MINUS_B = np.array([-0.385, -1.938, 0.995])
GOOD_A_B_AVG = np.array([3.9245, -0.834, -2.0205])
C_VEC = [24.117, -20.135, -52.518]
GOOD_A_MINUS_C = np.array([3.865, -5.918, 2.495])
GOOD_A_C_AVG = np.array([1.7995, 1.156, -2.7705])

VEC_1 = np.array([3.712, -1.585, -3.116])
VEC_2 = np.array([4.8760, -1.129, -3.265])
VEC_3 = np.array([5.498, -0.566, -2.286])
VEC_4 = np.array([5.464, -1.007, -0.948])

VEC_21 = np.array([-1.164, -0.456, 0.149])
VEC_23 = np.array([0.622, 0.563, 0.979])
VEC_34 = np.array([-0.034, -0.441, 1.338])

UNIT_VEC_3 = np.array([0.91922121129527656, -0.094630630337054641, -0.38220074372881085])
ANGLE_123 = 120.952786591
DIH_1234 = 39.4905248514

# -Radial Correction- #

CORR_KEY = 'corr'
COORD_KEY = 'coord'
FREE_KEY = 'free_energy'
RAD_KEY_SEQ = [COORD_KEY, FREE_KEY, CORR_KEY]


def expected_dir_data():
    """
    :return: The data structure that's expected from `find_files_by_dir`
    """
    return {os.path.abspath(os.path.join(FES_DIR, "1.00")): ['fes.out'],
            os.path.abspath(os.path.join(FES_DIR, "2.75")): ['fes.out', 'fes_cont.out'],
            os.path.abspath(os.path.join(FES_DIR, "5.50")): ['fes.out', 'fes_cont.out'],
            os.path.abspath(os.path.join(FES_DIR, "multi")): ['fes.out', 'fes_cont.out',
                                                              'fes_cont2.out', 'fes_cont3.out'],
            os.path.abspath(os.path.join(FES_DIR, "no_overwrite")): ['fes.out'], }


def csv_data():
    """
    :return: Test data as a list of dicts.
    """
    rows = [{CORR_KEY: 123.42, COORD_KEY: "75", FREE_KEY: True},
            {CORR_KEY: 999.43, COORD_KEY: "yellow", FREE_KEY: False}]
    return rows


def is_one_of_type(val, types):
    """Returns whether the given value is one of the given types.

    :param val: The value to evaluate
    :param types: A sequence of types to check against.
    :return: Whether the given value is one of the given types.
    """
    result = False
    val_type = type(val)
    for tt in types:
        if val_type is tt:
            result = True
    return result


# Tests #

class TestRateCalc(unittest.TestCase):
    """
    Tests calculation of a rate coefficient by the Eyring equation.
    """
    def test_calc_k(self):
        temp = 900.0
        delta_g = 53.7306
        rate_coeff = calc_k(temp, delta_g)
        self.assertEqual(rate_coeff, 1.648326791137026)

    # def test_calc_k_real(self):
    #     temp = 300.0
    #     delta_g = 36
    #     rate_coeff = calc_k(temp, delta_g)
    #     rate_coeff2 = calc_k(temp, 12.3)
    #     print(rate_coeff2/rate_coeff)
    #     print("Rate coefficient in s^-1: {}".format(rate_coeff))
    #     print("Timescale in s: {}".format(1/rate_coeff))
    #     print("Timescale in min: {}".format(1/rate_coeff/60))
    #     print("Timescale in hours: {}".format(1/rate_coeff/60/60))
    #     print("Timescale in days: {}".format(1/rate_coeff/60/60/24))
    #     print("Timescale in months: {}".format(1/rate_coeff/60/60/24/30))
    #     print("Timescale in years: {}".format(1/rate_coeff/60/60/24/365.25))
    #
    #
    # def test_calc_k_real2(self):
    #     temp = 300.0
    #     delta_g = 12.3
    #     rate_coeff = calc_k(temp, delta_g)
    #     print("Rate coefficient in s^-1: {}".format(rate_coeff))
    #     print("Timescale in s: {}".format(1/rate_coeff))
    #     print("Timescale in ms: {}".format(1000/rate_coeff))
    #     print("Timescale in microseconds: {}".format(1000*1000/rate_coeff))


class TestFindFiles(unittest.TestCase):
    """
    Tests for the file finder.
    """
    def test_find(self):
        found = find_files_by_dir(FES_DIR, DEF_FILE_PAT)
        exp_data = expected_dir_data()
        self.assertEqual(len(exp_data), len(found))
        for key, files in exp_data.items():
            found_files = found.get(key)
            try:
                self.assertEqual(len(files), len(found_files))
            except AttributeError:
                self.assertEqual(files, found_files)


class TestCheckFileFileList(unittest.TestCase):
    """
    Tests for the file finder.
    """
    def test_NoneOnly(self):
        try:
            found_list = check_for_files(None, None)
            self.assertFalse(found_list)
        except InvalidDataError as e:
            self.assertTrue("No files to process" in e.args[0])

    def test_NoSuchFile(self):
        try:
            found_list = check_for_files("ghost.com", None)
            self.assertFalse(found_list)
        except IOError as e:
            self.assertTrue("ghost.com" in e.args[0])

    def test_WrongInputType(self):
        try:
            found_list = check_for_files(None, 123)
            self.assertFalse(found_list)
        except InvalidDataError as e:
            self.assertTrue("file" in e.args[0])

    def test_name_only(self):
        found_list = check_for_files(ELEM_DICT_FILE, None)
        self.assertTrue(len(found_list) == 1)
        self.assertTrue(os.path.relpath(ELEM_DICT_FILE) == found_list[0])

    def testListFile(self):
        found_list = check_for_files(None, FILE_LIST)
        self.assertTrue(len(found_list) == 4)

    def testList(self):
        file_list = file_rows_to_list(FILE_LIST)
        found_list = check_for_files(None, file_list)
        self.assertTrue(len(found_list) == 4)

    def testListWithMissingFile(self):
        # found_list = check_for_files(None, FILE_LIST_W_MISSING_FILE)
        try:
            found_list = check_for_files(None, FILE_LIST_W_MISSING_FILE)
            self.assertFalse(found_list)
        except IOError as e:
            self.assertTrue("ghost.csv" in e.args[0])

    def testSearchCurrentDir(self):
        # this test assumes only only license file
        found_list = check_for_files(None, None, search_pattern="LICENSE")
        self.assertEqual(len(found_list), 1)
        self.assertEqual(os.path.relpath(found_list[0]), 'LICENSE')

    def testInvalidSearchDir(self):
        try:
            found_list = check_for_files(None, None, search_pattern="py", search_dir="ghost")
            self.assertFalse(found_list)
        except InvalidDataError as e:
            self.assertTrue("Could not find" in e.args[0])

    def testSearchSubDir(self):
        try:
            silent_remove(NEW_DIR, dir_with_files=True)
            make_dir(NEW_DIR)
            created_files = [TEST_DIR_SEARCH_FILE, TEST_DIR_SEARCH_FILE2]
            for fname in created_files:
                with open(fname, 'w') as f:
                    f.write("file is opened for business")
            found_list = check_for_files(None, None, search_pattern="unique", search_dir=SUB_DATA_DIR,
                                         search_sub_dir=True)
            self.assertEqual(found_list, [os.path.relpath(TEST_DIR_SEARCH_FILE2)])
        finally:
            silent_remove(NEW_DIR, disable=DISABLE_REMOVE, dir_with_files=True)

    def testSearchSubDirAltPat(self):
        try:
            silent_remove(NEW_DIR, dir_with_files=True)
            make_dir(NEW_DIR)
            created_files = [TEST_DIR_SEARCH_FILE, TEST_DIR_SEARCH_FILE2]
            for fname in created_files:
                with open(fname, 'w') as f:
                    f.write("file is opened for business")
            found_list = check_for_files(None, None, search_pattern="*unique", search_dir=SUB_DATA_DIR,
                                         search_sub_dir=True)
            self.assertEqual(found_list, [os.path.relpath(TEST_DIR_SEARCH_FILE2)])
        finally:
            silent_remove(NEW_DIR, disable=DISABLE_REMOVE, dir_with_files=True)

    def testSearchButDoNotFindDir(self):
        try:
            silent_remove(NEW_DIR, dir_with_files=True)
            found_list = check_for_files(None, None, search_pattern="unique", search_dir=SUB_DATA_DIR)
            self.assertFalse(found_list)
        except InvalidDataError as e:
            self.assertTrue("No files to process" in e.args[0])
        finally:
            silent_remove(NEW_DIR, disable=DISABLE_REMOVE, dir_with_files=True)

    def testSearchButDoNotFindSubDir(self):
        try:
            silent_remove(NEW_DIR, dir_with_files=True)
            found_list = check_for_files(None, None, search_pattern="unique", search_dir=SUB_DATA_DIR,
                                         search_sub_dir=True)
            self.assertFalse(found_list)
        except InvalidDataError as e:
            self.assertTrue("No files to process" in e.args[0])
        finally:
            silent_remove(NEW_DIR, disable=DISABLE_REMOVE, dir_with_files=True)


class TestMakeDir(unittest.TestCase):
    def testExistingDir(self):
        try:
            hello = make_dir(SUB_DATA_DIR)
            self.assertTrue(hello is None)
        except NotFoundError:
            self.fail("make_dir() raised NotFoundError unexpectedly!")

    def testNewDir(self):
        try:
            silent_remove(NEW_DIR)
            make_dir(NEW_DIR)
            self.assertTrue(os.path.isdir(NEW_DIR))
        finally:
            silent_remove(NEW_DIR, disable=DISABLE_REMOVE)


class TestReadJson(unittest.TestCase):
    def testJsonToDict(self):
        good_project_dict = {'EXPECTED_BASIS': 'def2tzvp',
                             'LOG_DIR': 'tests/test_data/gaussian_output_condensed_phase',
                             'RXN_STEPS': {'1': {'TYPES': ['ESTER_H_OVERALL'], 'NAMES': ['TPA', 'IPA'],
                                                 'STEPS': ['INI_REACTS', 'STEP1_REACTS', 'STEP1_TS', 'STEP1_PRODS']}},
                             'SINGLE_TEMP_K': 433.15, 'TEMP_RANGE': '353.15,463.15,10', 'TEMPS_K': [363.15, 423.15],
                             'VIB_SCALE': 0.9505}
        j_file = os.path.join(SUB_DATA_DIR, "project_info.json")
        project_dict = read_json(j_file)
        self.assertEqual(project_dict, good_project_dict)

    def testJsonToDictCatchError(self):
        try:
            j_file = os.path.join(SUB_DATA_DIR, "project_info_error.json")
            read_json(j_file)
            self.assertFalse("I should not be reached")
        except InvalidDataError as e:
            self.assertTrue("Error in reading JSON format" in e.args[0])


class TestReadPDB(unittest.TestCase):
    def testReadPDBFname(self):
        pdb_fname = os.path.join(SUB_DATA_DIR, "pet.pdb")
        pdb_data_dict = process_pdb_file(pdb_fname)
        self.assertEqual(len(pdb_data_dict), 4)
        self.assertEqual(pdb_data_dict[NUM_ATOMS], 33)
        self.assertEqual(pdb_data_dict[SEC_HEAD], ['COMPND    PET FRAGMENT 2',
                                                   'AUTHOR    GENERATED BY OPEN BABEL 2.3.90'])
        self.assertEqual(pdb_data_dict[SEC_ATOMS][-1], ['HETATM', '   33', ' O    ', 'UNL  ', 1,
                                                        -3.539, 1.04, 10.096, '  1.00  0.00           O'])
        self.assertEqual(len(pdb_data_dict[SEC_TAIL]), 35)
        self.assertEqual(pdb_data_dict[SEC_TAIL][-3], 'CONECT   33   26   29')

    def testReadPDBFnameAtomInfoOnly(self):
        pdb_fname = os.path.join(SUB_DATA_DIR, "pet.pdb")
        pdb_data_dict = process_pdb_file(pdb_fname, atom_info_only=True)
        self.assertEqual(len(pdb_data_dict), 4)
        self.assertEqual(pdb_data_dict[NUM_ATOMS], 33)
        self.assertEqual(pdb_data_dict[SEC_HEAD], ['COMPND    PET FRAGMENT 2',
                                                   'AUTHOR    GENERATED BY OPEN BABEL 2.3.90'])
        self.assertEqual(pdb_data_dict[SEC_ATOMS][33][ATOM_TYPE], ' O')
        self.assertTrue(np.allclose(pdb_data_dict[SEC_ATOMS][33][ATOM_COORDS], np.array([-3.539, 1.04, 10.096])))
        self.assertEqual(len(pdb_data_dict[SEC_TAIL]), 35)
        self.assertEqual(pdb_data_dict[SEC_TAIL][-3], 'CONECT   33   26   29')


class TestReadFirstRow(unittest.TestCase):
    def testFirstRow(self):
        self.assertListEqual(CSV_HEADER, read_csv_header(CSV_FILE))

    def testSecondRow(self):
        second_row = read_csv_header(ALT_CSV_FILE, return_second_row=True)
        self.assertEqual(second_row, read_csv_header(CSV_FILE))

    def testEmptyFile(self):
        self.assertIsNone(read_csv_header(EMPTY_CSV))


class TestIOMethods(unittest.TestCase):
    def testReadTplNoSuchTpl(self):
        try:
            tpl_str = read_tpl("ghost.tpl")
            self.assertFalse(tpl_str)
        except TemplateNotReadableError as e:
            self.assertTrue("Couldn't read template at: 'ghost.tpl'" in e.args[0])

    def testMakeDir(self):
        # provide a file name not a dir
        try:
            make_dir(ELEM_DICT_FILE)
            # should raise exception before next line
            self.assertFalse(True)
        except NotFoundError as e:
            self.assertTrue("Resource exists and is not a dir" in e.args[0])

    def testFileRowsToList(self):
        # this function should skip blank lines
        test_rows = file_rows_to_list(FILE_LIST)
        good_rows = ['tests/test_data/common/diff_lines_base_file.csv',
                     'tests/test_data/common/diff_lines_miss_line.csv',
                     'tests/test_data/common/diff_lines_miss_val.csv',
                     'tests/test_data/common/diff_lines_one_nan.csv']
        self.assertTrue(test_rows == good_rows)


class TestArgParseFunctions(unittest.TestCase):
    def testOverwriteConfigVals(self):
        arg_config_keyword_dict = {'file': 'file',
                                   'list_fname': 'list_fname',
                                   'out_dir': 'out_dir',
                                   'float_val': 'float_val',
                                   'new_fname': 'new_fname',
                                   'add_elements': 'add_elements'}
        default_val_dict = {'list_fname': None,
                            'file': None,
                            'out_dir': None,
                            'float_val': 2.3,
                            'new_fname': None,
                            'pdb_print_format': '{:6s}{:>5}{:^6s}{:5s}{:>4}    {:8.3f}{:8.3f}{:8.3f}{:22s}{:>2s}{:s}',
                            'add_elements': False,
                            }

        args = Namespace(add_elements=True,
                         file="do_not_actually_exist.txt",
                         float_val=9.8,
                         config=default_val_dict.copy(),
                         list_fname=None,
                         new_fname=None,
                         out_dir=None)
        args.config['new_fname'] = "hello_world.txt"

        overwrite_config_vals(args, arg_config_keyword_dict, default_val_dict)
        # check that did not overwrite config value with default arg val
        self.assertEqual(args.config['new_fname'], "hello_world.txt")
        # check that did overwrite add_elements flag
        self.assertTrue(args.config['add_elements'])
        # check float comparision
        self.assertAlmostEqual(9.8, args.config['float_val'])


class TestRoundingToSigFig(unittest.TestCase):
    def testRoundTo12thDecimal(self):
        # helps in printing, so files aren't different only due to expected machine precision (:.12f, but keep as float)
        result = round_to_12th_decimal(8.76541113456789012345)
        good_result = 8.765411134568
        self.assertTrue(result == good_result)

    def testStandardUse(self):
        self.assertAlmostEqual(round_sig_figs(1111.111111111), 1111.11)

    def testSpecifySigFigs(self):
        self.assertAlmostEqual(round_sig_figs(1111.111111111, sig_figs=3), 1110.0)

    def testIntSpecifySigFigs(self):
        rounded_val = round_sig_figs(111111, sig_figs=3)
        self.assertIsInstance(rounded_val, int)
        self.assertEqual(rounded_val, 111000)

    def testFloat64RoundSigFigs(self):
        input_num = np.float64(1111.111111111)
        out_num = round_sig_figs(input_num, sig_figs=3)
        self.assertEqual(type(input_num), type(out_num))

    def testFloat32RoundSigFigs(self):
        input_num = np.float32(1111.111111111)
        out_num = round_sig_figs(input_num, sig_figs=3)
        self.assertEqual(type(input_num), type(out_num))


class TestRoundToFraction(unittest.TestCase):
    def testClosestInteger(self):
        test_array = [0.256, 0.02, 1.8767, 2.432, 8.174, 4.63, 2.82, 5.311]
        expected_out = [0., 0., 2., 2., 8., 5., 3., 5.]
        array_out = round_to_fraction(test_array, 1)
        self.assertTrue(np.allclose(array_out, expected_out))

    def testClosestIntegerFloat(self):
        test_array = [0.256, 0.02, 1.8767, 2.432, 8.174, 4.63, 2.82, 5.311]
        expected_out = [0., 0., 2., 2., 8., 5., 3., 5.]
        array_out = round_to_fraction(test_array, 1.0)
        self.assertTrue(np.allclose(array_out, expected_out))

    def testClosestTenth(self):
        test_array = [0.256, 0.02, 1.8767, 2.432, 8.174, 4.63, 2.82, 5.311]
        expected_out = [0.3, 0.0, 1.9, 2.4, 8.2, 4.6, 2.8, 5.3]
        array_out = round_to_fraction(test_array, 0.1)
        self.assertTrue(np.allclose(array_out, expected_out))

    def testClosestHundredth(self):
        test_array = [0.25612, 9., 1.8767, 2.4323145, 8.174048, 4.638, 2.82, 5.31111]
        expected_out = [0.26, 9., 1.88, 2.43, 8.17, 4.64, 2.82, 5.31]
        array_out = round_to_fraction(test_array, 0.01)
        self.assertTrue(np.allclose(array_out, expected_out))

    def testClosestMilli(self):
        test_array = [0.25612, 9., 1.8767, 2.4323145, 8.174048, 4.638, 2.82, 5.31111]
        expected_out = [0.256, 9., 1.877, 2.432, 8.174, 4.638, 2.82, 5.311]
        array_out = round_to_fraction(test_array, 0.001)
        self.assertTrue(np.allclose(array_out, expected_out))

    def testClosestHalf(self):
        test_array = [0.256, 0.02, 1.8767, 2.432, 8.174, 4.63, 2.82, 5.311]
        expected_out = [0.5, 0.0, 2.0, 2.5, 8.0, 4.5, 3.0, 5.5]
        array_out = round_to_fraction(test_array, 0.5)
        self.assertTrue(np.allclose(array_out, expected_out))

    def testClosestTwentieth(self):
        test_array = [0.256, 0.02, 1.8767, 2.432, 8.174, 4.63, 2.82, 5.311]
        expected_out = [0.25, 0., 1.9, 2.45, 8.15, 4.65, 2.8, 5.3]
        array_out = round_to_fraction(test_array, 0.05)
        self.assertTrue(np.allclose(array_out, expected_out))

    def testClosestQuarter(self):
        test_array = [0.256, 0.02, 1.8767, 2.432, 8.174, 4.63, 2.82, 5.311]
        expected_out = [0.25, 0.00, 2.00, 2.50, 8.25, 4.75, 2.75, 5.25]
        array_out = round_to_fraction(test_array, 0.25)
        self.assertTrue(np.allclose(array_out, expected_out))

    def testClosestThird(self):
        test_array = [0.256, 0.02, 1.8767, 2.432, 8.174, 4.63, 2.82, 5.311]
        expected_out = [0.33, 0.00, 2.00, 2.33, 8.33, 4.67, 2.67, 5.33]
        array_out = round_to_fraction(test_array, 1./3.)
        array_out = np.around(array_out, 2)
        self.assertTrue(np.allclose(array_out, expected_out))

    def testNotFractionOfOne(self):
        try:
            test_array = [0.256, 0.02, 1.8767, 2.432, 8.174, 4.63, 2.82, 5.311]
            round_to_fraction(test_array, 0.0011)
            self.assertFalse("I should have had an error before I got here.")
        except InvalidDataError as e:
            self.assertTrue("evenly divides into 1" in e.args[0])

    def testFractionTen(self):
        try:
            test_array = [0.256, 0.02, 1.8767, 2.432, 8.174, 4.63, 2.82, 5.311]
            round_to_fraction(test_array, 10.)
            self.assertFalse("I should have had an error before I got here.")
        except InvalidDataError as e:
            self.assertTrue("evenly divides into 1" in e.args[0])

    def testFractionGreaterThanOne(self):
        try:
            test_array = [0.256, 0.02, 1.8767, 2.432, 8.174, 4.63, 2.82, 5.311]
            round_to_fraction(test_array, 11.1)
            self.assertFalse("I should have had an error before I got here.")
        except InvalidDataError as e:
            self.assertTrue("evenly divides into 1" in e.args[0])

    def testNegNumber(self):
        try:
            test_array = [0.256, 0.02, 1.8767, 2.432, 8.174, 4.63, 2.82, 5.311]
            round_to_fraction(test_array, -.1)
            self.assertFalse("I should have had an error before I got here.")
        except InvalidDataError as e:
            self.assertTrue("positive number" in e.args[0])

    def testNonFloat(self):
        try:
            test_array = [0.256, 0.02, 1.8767, 2.432, 8.174, 4.63, 2.82, 5.311]
            round_to_fraction(test_array, "0.1")
            self.assertFalse("I should have had an error before I got here.")
        except InvalidDataError as e:
            self.assertTrue("positive number" in e.args[0])


class TestNaturalSorting(unittest.TestCase):
    def testTextFloats(self):
        float_str_list = ["12.1", "94.8", "6.2", "36.1"]
        float_str_list.sort(key=natural_keys)
        good_sorted_str_list = ['6.2', '12.1', '36.1', '94.8']
        self.assertEqual(float_str_list, good_sorted_str_list)

    def testMolFormulas(self):
        """
        when I had CH4, I added a trick: after this sort, sort by length:
            dict_keys.sort(key=natural_keys)
            dict_keys.sort(key=len)
        that was good enough for me!
        """
        str_list = ["C2H6", "C6H6", "C10H20", "C2H4"]
        str_list.sort(key=natural_keys)
        good_sorted_str_list = ["C2H4", "C2H6", "C6H6", "C10H20"]
        self.assertEqual(str_list, good_sorted_str_list)

    def testTextAlphaNumeric(self):
        str_list = ["2018AnotherFile.txt", "1956OlderFile.txt", "2019-02-11File_name.txt", "422MuchOlder.text"]
        str_list.sort(key=natural_keys)
        sorted_str_list = ["422MuchOlder.text", "1956OlderFile.txt", "2018AnotherFile.txt", "2019-02-11File_name.txt"]
        self.assertEqual(str_list, sorted_str_list)


class TestUniqueList(unittest.TestCase):
    def testSameInOut(self):
        unique_sorted_list = ["alligator", "elephant", "walrus"]
        returned_list = unique_list(unique_sorted_list)
        self.assertEqual(returned_list, unique_sorted_list)

    def testRepeatIn(self):
        unique_sorted_list = ["alligator", "elephant", "walrus"]
        repeat_sorted_list = ["alligator", "elephant", "elephant", "elephant", "walrus", "walrus"]
        returned_list = unique_list(repeat_sorted_list)
        self.assertEqual(returned_list, unique_sorted_list)

    def testUnsortedTuple(self):
        list_unique = ["walrus", "alligator", "elephant"]
        repeat_in_tuple = ("walrus", "alligator", "walrus", "elephant", "walrus")
        returned_list = unique_list(repeat_in_tuple)
        self.assertEqual(returned_list, list_unique)


class TestFnameManipulation(unittest.TestCase):
    def testOutFname(self):
        """
        Check for prefix addition.
        """
        self.assertTrue(create_out_fname(ORIG_WHAM_PATH, prefix=OUT_PFX).endswith(
            os.sep + OUT_PFX + ORIG_WHAM_FNAME))

    def testOutFnameRemovePrefix(self):
        """
        Check for prefix addition after prefix removal.
        """
        prefix_to_remove = 'ghost'
        beginning_name = os.path.join(DATA_DIR, prefix_to_remove + ORIG_WHAM_FNAME)
        good_end_name = os.path.join(DATA_DIR, OUT_PFX + ORIG_WHAM_FNAME)
        new_name = create_out_fname(beginning_name, prefix=OUT_PFX, remove_prefix=prefix_to_remove)
        self.assertTrue(new_name == good_end_name)

    def testOutFnameGiveExt(self):
        beginning_name = os.path.join(DATA_DIR, ORIG_WHAM_ROOT)
        good_end_name = os.path.join(DATA_DIR, ORIG_WHAM_FNAME)
        new_name = create_out_fname(beginning_name, ext='.txt')
        self.assertTrue(new_name == good_end_name)

    def testOutFnameGiveExtNoPeriod(self):
        beginning_name = os.path.join(DATA_DIR, ORIG_WHAM_ROOT)
        good_end_name = os.path.join(DATA_DIR, ORIG_WHAM_FNAME)
        new_name = create_out_fname(beginning_name, ext='txt')
        self.assertTrue(new_name == good_end_name)

    def testOutFnameNoPath(self):
        """
        Check for prefix addition, without a path in the input
        """
        returned_name = create_out_fname(ORIG_WHAM_ROOT, prefix=OUT_PFX, ext="txt")
        self.assertTrue(returned_name.endswith(os.sep + OUT_PFX + ORIG_WHAM_FNAME))

    def testGetRootName(self):
        root_name = get_fname_root(ORIG_WHAM_PATH)
        self.assertEqual(root_name, ORIG_WHAM_ROOT)
        self.assertNotEqual(root_name, ORIG_WHAM_FNAME)
        self.assertNotEqual(root_name, ORIG_WHAM_PATH)

    def testRelPath(self):
        input_fname = os.path.join(DATA_DIR, ORIG_WHAM_ROOT)
        good_returned_name = os.path.relpath(os.path.join(DATA_DIR, ORIG_WHAM_FNAME))
        returned_name = create_out_fname(input_fname, ext='txt', rel_path=True)
        self.assertTrue(returned_name == good_returned_name)


class TestReadCsvDict(unittest.TestCase):
    def testReadAtomNumDict(self):
        # Will renumber atoms and then sort them
        test_dict = read_csv_dict(ATOM_DICT_FILE)
        self.assertEqual(test_dict, GOOD_ATOM_DICT)

    def testReadAtomNumDictAsFloat(self):
        # Will renumber atoms and then sort them
        atom_float_dict = {}
        for key, val in GOOD_ATOM_DICT.items():
            atom_float_dict[str(key)] = float(val)
        test_dict = read_csv_dict(ATOM_DICT_FILE, str_float=True)
        self.assertEqual(test_dict, atom_float_dict)

    def testReadPDBDict(self):
        test_type = ' HY1'
        test_elem = ' H'
        test_dict = read_csv_dict(ELEM_DICT_FILE, pdb_dict=True)
        self.assertTrue(test_type in test_dict)
        self.assertEqual(test_elem, test_dict[test_type])
        self.assertEqual(31, len(test_dict))

    def testReadPDBDictWithBlanks(self):
        element_dict_file = os.path.join(SUB_DATA_DIR, "element_dict_w_blanks.csv")
        test_type = 'HNT1'
        test_elem = ' H'
        test_dict = read_csv_dict(element_dict_file, pdb_dict=True)
        self.assertTrue(test_type in test_dict)
        self.assertEqual(test_elem, test_dict[test_type])
        self.assertEqual(31, len(test_dict))

    def testStringDictAsInt(self):
        # Check that fails elegantly by passing returning value error
        try:
            test_dict = read_csv_dict(ELEM_DICT_FILE, one_to_one=False)
            self.assertFalse(test_dict)
        except ValueError as e:
            self.assertTrue("invalid literal for int()" in e.args[0])

    def testStringDictCheckDups(self):
        # Check that fails elegantly
        try:
            test_dict = read_csv_dict(ELEM_DICT_FILE, ints=False, )
            self.assertFalse(test_dict)
        except InvalidDataError as e:
            self.assertTrue("Did not find a 1:1 mapping" in e.args[0])

    def testStringDictExtraTerm(self):
        # Check that fails elegantly
        input_dict = os.path.join(SUB_DATA_DIR, "element_dict_extra_term.csv")
        try:
            test_dict = read_csv_dict(input_dict, ints=False, )
            self.assertFalse(test_dict)
        except InvalidDataError as e:
            self.assertTrue("exactly two" in e.args[0])

    def testStringDictDuplicateTerm(self):
        # Check that fails elegantly
        input_dict = os.path.join(SUB_DATA_DIR, "element_dict_dup.csv")
        try:
            test_dict = read_csv_dict(input_dict, ints=False, )
            self.assertFalse(test_dict)
        except InvalidDataError as e:
            self.assertTrue("non-unique" in e.args[0])

    def testElementDictKeyTooLong(self):
        # Check that fails elegantly
        input_dict = os.path.join(SUB_DATA_DIR, "element_dict_too_long_key.csv")
        try:
            test_dict = read_csv_dict(input_dict, pdb_dict=True)
            self.assertFalse(test_dict)
        except InvalidDataError as e:
            self.assertTrue("no more than" in e.args[0])


class TestReadCsv(unittest.TestCase):
    def testReadCsv(self):
        """
        Verifies the contents of the CSV file.
        """
        result = read_csv(CSV_FILE)
        self.assertTrue(result)
        for row in result:
            self.assertEqual(3, len(row))
            self.assertIsNotNone(row.get(FREE_KEY, None))
            self.assertIsInstance(row[FREE_KEY], str)
            self.assertIsNotNone(row.get(CORR_KEY, None))
            self.assertIsInstance(row[CORR_KEY], str)
            self.assertIsNotNone(row.get(COORD_KEY, None))
            self.assertIsInstance(row[COORD_KEY], str)

    def testReadTypedCsvAllConv(self):
        """
        Verifies the contents of the CSV file using the all_conv function.
        """
        result = read_csv(CSV_FILE, all_conv=float)
        self.assertTrue(result)
        for row in result:
            self.assertEqual(3, len(row))
            self.assertIsNotNone(row.get(FREE_KEY, None))
            self.assertTrue(is_one_of_type(row[FREE_KEY], FRENG_TYPES))
            self.assertIsNotNone(row.get(CORR_KEY, None))
            self.assertTrue(is_one_of_type(row[CORR_KEY], FRENG_TYPES))
            self.assertIsNotNone(row.get(COORD_KEY, None))
            self.assertIsInstance(row[COORD_KEY], float)

    def testReadCsvToList(self):
        read_list, header = read_csv_to_list(BOX_SIZES_HEADER_COMMA_SEP, header=True)
        self.assertEqual(header, ['x', 'y', 'z'])
        self.assertTrue(len(read_list), 4)
        self.assertTrue(len(read_list[-1]), 3)

    def testReadCsvToListNoHeader(self):
        read_list, header = read_csv_to_list(BOX_SIZES_NO_HEADER_COMMA_SEP)
        self.assertEqual(header, [])
        self.assertTrue(len(read_list), 4)
        self.assertTrue(len(read_list[-1]), 3)


class TestWriteCsv(unittest.TestCase):
    def testWriteCsv(self):
        tmp_dir = None
        data = csv_data()
        try:
            tmp_dir = tempfile.mkdtemp()
            tgt_fname = create_out_fname(SHORT_WHAM_PATH, prefix=OUT_PFX, base_dir=tmp_dir)
            # write_csv(data, tgt_fname, RAD_KEY_SEQ)
            with capture_stdout(write_csv, data, tgt_fname, RAD_KEY_SEQ) as output:
                self.assertTrue("Wrote file:" in output)
            csv_result = read_csv(tgt_fname,
                                  data_conv={FREE_KEY: str_to_bool,
                                             CORR_KEY: float,
                                             COORD_KEY: str, })
            self.assertEqual(len(data), len(csv_result))
            for i, csv_row in enumerate(csv_result):
                self.assertDictEqual(data[i], csv_row)
        finally:
            shutil.rmtree(tmp_dir)

    def testAppendCsv(self):
        tmp_dir = None
        data = csv_data()
        try:
            tmp_dir = tempfile.mkdtemp()
            tgt_fname = create_out_fname(SHORT_WHAM_PATH, prefix=OUT_PFX, base_dir=tmp_dir)
            # write_csv(data, tgt_fname, RAD_KEY_SEQ)
            with capture_stdout(write_csv, data, tgt_fname, RAD_KEY_SEQ, mode="a") as output:
                self.assertTrue("Appended:" in output)
            csv_result = read_csv(tgt_fname,
                                  data_conv={FREE_KEY: str_to_bool,
                                             CORR_KEY: float,
                                             COORD_KEY: str, })
            dict_from_reading_append = [{str(data[0][FREE_KEY]): str(data[1][FREE_KEY]),
                                        str(data[0][CORR_KEY]): str(data[1][CORR_KEY]),
                                        data[0][COORD_KEY]: data[1][COORD_KEY], }]
            self.assertEqual(len(dict_from_reading_append), len(csv_result))
            for i, csv_row in enumerate(csv_result):
                self.assertDictEqual(dict_from_reading_append[i], csv_row)
        finally:
            shutil.rmtree(tmp_dir)

    def testRoundNum(self):
        # like testWriteCsv, but have it round away extra digits
        tmp_dir = None
        data = csv_data()
        for data_dict in data:
            data_dict[CORR_KEY] += 0.0024
        try:
            tmp_dir = tempfile.mkdtemp()
            tgt_fname = create_out_fname(SHORT_WHAM_PATH, prefix=OUT_PFX, base_dir=tmp_dir)
            write_csv(data, tgt_fname, RAD_KEY_SEQ, round_digits=2)
            csv_result = read_csv(tgt_fname,
                                  data_conv={FREE_KEY: str_to_bool,
                                             CORR_KEY: float,
                                             COORD_KEY: str, })
            self.assertEqual(len(data), len(csv_result))
            for i, csv_row in enumerate(csv_result):
                self.assertDictEqual(data[i], csv_row)
        finally:
            shutil.rmtree(tmp_dir)

    def testWriteCsvToStdOut(self):
        data = csv_data()
        good_output_list = ['"coord","free_energy","corr"',
                            '"75",True,123.42',
                            '"yellow",False,999.43', '']
        with capture_stdout(print_csv_stdout, data, RAD_KEY_SEQ) as output:
            output_list = output.split('\r\n')
            self.assertTrue(output_list == good_output_list)


class TestNDARRAYFromFile(unittest.TestCase):
    def testReadSpaceSeparated(self):
        test_array = np_float_array_from_file(BOX_SIZES_NO_HEADER_SPACE_SEP)
        self.assertTrue(np.allclose(test_array, GOOD_BOX_NDARRAY))

    def testReadSpaceSeparatedWithHeaderFlagged(self):
        good_header_row = ['#', 'x', 'y', 'z']
        test_array, test_header = np_float_array_from_file(BOX_SIZES_HEADER_SPACE_SEP, header=True)
        self.assertEqual(test_header, good_header_row)
        self.assertTrue(np.allclose(test_array, GOOD_BOX_NDARRAY))

    def testReadCommaSeparated(self):
        test_array = np_float_array_from_file(BOX_SIZES_NO_HEADER_COMMA_SEP, delimiter=',')
        self.assertTrue(np.allclose(test_array, GOOD_BOX_NDARRAY))

    def testReadCommaSeparatedWithHeaderNotFlagged(self):
        with capture_stderr(np_float_array_from_file, BOX_SIZES_HEADER_COMMA_SEP, delimiter=',') as output:
            self.assertTrue("header" in output)
        test_array = np_float_array_from_file(BOX_SIZES_HEADER_COMMA_SEP, delimiter=',')
        self.assertTrue(np.allclose(test_array, GOOD_BOX_NDARRAY_ROW_NAN, equal_nan=True))

    def testReadCommaSeparatedWithHeaderFlagged(self):
        good_header_row = ['x', 'y', 'z']
        test_array, test_header = np_float_array_from_file(BOX_SIZES_HEADER_COMMA_SEP, delimiter=',', header=True)
        self.assertEqual(test_header, good_header_row)
        self.assertTrue(np.allclose(test_array, GOOD_BOX_NDARRAY))

    def testReadNDArrayVectorError(self):
        caught_error = False
        vector_vals = os.path.join(SUB_DATA_DIR, 'vector_input.txt')
        try:
            np_float_array_from_file(vector_vals)
        except InvalidDataError as e:
            caught_error = True
            self.assertTrue('File contains a vector' in e.args[0])
        self.assertTrue(caught_error)

    def testReadNDArrayValueError(self):
        # TODO: copy this but test and use hist option, to improve coverage
        with capture_stderr(np_float_array_from_file, FLOAT_AND_NON, header=1, delimiter=",") as output:
            self.assertTrue("'nan' will be returned" in output)
        good_header = ['pka_203', '(0, 1)', '(0, 1)_max_rls', '(0, -1)_max_rls', '(0, 1)_max_path',
                       '(0, 1)_max_path_flow']
        good_data_array = np.asarray([[6.10918105, 1.04301557, np.nan, np.nan, np.nan, 0.91252645],
                                      [4.33909619, 1.09081880, np.nan, np.nan, np.nan, 0.87450673],
                                      [5.54534891, 1.06042369, np.nan, np.nan, np.nan, 0.65756597],
                                      [5.29792317, 1.08224906, np.nan, np.nan, np.nan, 0.80011857],
                                      [5.99200576, 1.06021529, np.nan, 0., np.nan, 0.48421892]])
        data_array, header_row = np_float_array_from_file(FLOAT_AND_NON, header=1, delimiter=",")
        self.assertTrue(header_row == good_header)
        self.assertTrue(np.allclose(data_array, good_data_array, equal_nan=True))

    def testReadNDArrayValueErrorNoHeaderSpecified(self):
        with capture_stderr(np_float_array_from_file, FLOAT_AND_NON, delimiter=",") as output:
            self.assertTrue("'nan' will be returned" in output)
        good_header = ['pka_203', '(0, 1)', '(0, 1)_max_rls', '(0, -1)_max_rls', '(0, 1)_max_path',
                       '(0, 1)_max_path_flow']
        good_data_array = np.asarray([[6.10918105, 1.04301557, np.nan, np.nan, np.nan, 0.91252645],
                                      [4.33909619, 1.09081880, np.nan, np.nan, np.nan, 0.87450673],
                                      [5.54534891, 1.06042369, np.nan, np.nan, np.nan, 0.65756597],
                                      [5.29792317, 1.08224906, np.nan, np.nan, np.nan, 0.80011857],
                                      [5.99200576, 1.06021529, np.nan, 0., np.nan, 0.48421892]])
        data_array, header_row = np_float_array_from_file(FLOAT_AND_NON, header=1, delimiter=",")
        self.assertTrue(header_row == good_header)
        self.assertTrue(np.allclose(data_array, good_data_array, equal_nan=True))


class TestListToFile(unittest.TestCase):
    def testWriteFormatListAppendList(self):
        try:
            list_of_lists_of_floats = [[1.2, 4.666698], [7.098, 89.1275]]
            list_to_file(list_of_lists_of_floats, LIST_OUT, list_format="{:6.2f} {:6.2f}")
            self.assertFalse(diff_lines(LIST_OUT, GOOD_FORMAT_LIST_OUT))
        finally:
            silent_remove(LIST_OUT, disable=DISABLE_REMOVE)
            pass

    def testWriteAppendList(self):
        list_of_strings = ['hello', 'friends']
        list_of_lists = [VEC_23, VEC_34]

        try:
            # list_to_file(list_of_strings, LIST_OUT)
            with capture_stdout(list_to_file, list_of_strings, LIST_OUT) as output:
                self.assertTrue("Wrote file: tests/test_data/common/temp.txt" in output)
            # list_to_file(VEC_21, LIST_OUT, mode="a")
            with capture_stdout(list_to_file, VEC_21, LIST_OUT, mode="a") as output:
                self.assertTrue("  Appended: tests/test_data/common/temp.txt" in output)
            # list_to_file(list_of_strings, LIST_OUT)
            with capture_stdout(list_to_file, list_of_lists, LIST_OUT, mode="a", print_message=False) as output:
                self.assertTrue(len(output) == 0)
            self.assertFalse(diff_lines(LIST_OUT, GOOD_LIST_OUT))
        finally:
            silent_remove(LIST_OUT, disable=DISABLE_REMOVE)
            pass

    def testWriteAppendStr(self):
        try:
            rel_list_out = os.path.relpath(LIST_OUT)
            string1 = "hello\nfriends\n"
            string2 = "-1.164\n-0.456\n0.149\n0.622 0.563 0.979\n-0.034 -0.441 1.338\n"
            # str_to_file(string1, rel_list_out, print_info=True)
            with capture_stdout(str_to_file, string1, rel_list_out, print_info=True) as output:
                self.assertTrue("Wrote file: tests/test_data/common/temp.txt" in output)
            # str_to_file(string2, rel_list_out, mode='a', print_info=True)
            with capture_stdout(str_to_file, string2, rel_list_out, mode="a", print_info=True) as output:
                self.assertTrue("  Appended: tests/test_data/common/temp.txt" in output)

            self.assertFalse(diff_lines(LIST_OUT, GOOD_LIST_OUT))
        finally:
            silent_remove(LIST_OUT, disable=DISABLE_REMOVE)
            pass

    def testListToCSV(self):
        try:
            good_out = os.path.join(SUB_DATA_DIR, "list_to_csv_good.txt")
            list_of_lists_of_floats = [[1.2, 4.666698], [7.098, 89.1275]]
            list_to_csv(list_of_lists_of_floats, LIST_OUT)
            self.assertFalse(diff_lines(good_out, LIST_OUT))
        finally:
            silent_remove(LIST_OUT, disable=DISABLE_REMOVE)
            pass

    def testListToCSVNotNested(self):
        try:
            good_out = os.path.join(SUB_DATA_DIR, "list_to_csv_not_nested_good.txt")
            list_of_floats = [1.2, 4.666698]
            list_to_csv(list_of_floats, LIST_OUT, round_digits=3)
            self.assertFalse(diff_lines(good_out, LIST_OUT))
        finally:
            silent_remove(LIST_OUT, disable=DISABLE_REMOVE)
            pass

    def testListToCSVRoundDigits(self):
        try:
            good_out = os.path.join(SUB_DATA_DIR, "list_to_csv_round_digits_good.txt")
            list_of_lists_of_floats = [[1.2, 4.666698], [7.098, 89.1275]]
            list_to_csv(list_of_lists_of_floats, LIST_OUT, round_digits=2)
            self.assertFalse(diff_lines(good_out, LIST_OUT))
        finally:
            silent_remove(LIST_OUT, disable=DISABLE_REMOVE)
            pass

    def testListToCSVListOfListsOfLists(self):
        list_of_lists_of_lists_of_floats = [[3.4], [[1.2, 4.666698], [7.098, 89.1275]]]
        try:
            good_out = os.path.join(SUB_DATA_DIR, "list_to_csv_extra_nesting_good.txt")
            list_to_csv(list_of_lists_of_lists_of_floats, LIST_OUT, round_digits=3)
            self.assertFalse(diff_lines(good_out, LIST_OUT))
        finally:
            silent_remove(LIST_OUT, disable=DISABLE_REMOVE)
            pass

    def testListToCSVGivenDict(self):
        input_dict = {0: [1.2, 4.666698]}
        try:
            list_to_csv(input_dict, LIST_OUT, round_digits=3)
            self.assertFalse("if you see me, raise and error!")
        except InvalidDataError as e:
            self.assertTrue("Expected" in e.args[0])
        finally:
            silent_remove(LIST_OUT, disable=DISABLE_REMOVE)
            pass


class TestFormatData(unittest.TestCase):
    def testFormatRows(self):
        raw = [{"a": 1.3333322333, "b": 999.222321}, {"a": 333.44422222, "b": 17.121}]
        fmt_std = [{'a': '1.3333', 'b': '999.2223'}, {'a': '333.4442', 'b': '17.1210'}]
        self.assertListEqual(fmt_std, fmt_row_data(raw, "{0:.4f}"))


class TestDiffLines(unittest.TestCase):
    def testSameFile(self):
        test_input = (DIFF_LINES_BASE_FILE, DIFF_LINES_BASE_FILE)
        diff_line_list = diff_lines(*test_input)
        self.assertEqual(len(diff_line_list), 0)
        with capture_stderr(diff_lines, *test_input) as output:
            self.assertFalse(output)

    def testMachinePrecDiff(self):
        test_input = (DIFF_LINES_BASE_FILE, DIFF_LINES_PREC_DIFF)
        diff_line_list = diff_lines(*test_input)
        self.assertEqual(len(diff_line_list), 0)
        with capture_stderr(diff_lines, *test_input) as output:
            self.assertTrue("floating point precision" in output)

    def testMachinePrecDiff2(self):
        test_input = (DIFF_LINES_PREC_DIFF, DIFF_LINES_BASE_FILE)
        diff_line_list = diff_lines(*test_input)
        self.assertEqual(len(diff_line_list), 0)
        with capture_stderr(diff_lines, *test_input) as output:
            self.assertTrue("floating point precision" in output)

    def testDiff(self):
        diffs = diff_lines(DIFF_LINES_ONE_VAL_DIFF, DIFF_LINES_BASE_FILE)
        self.assertEqual(len(diffs), 2)

    def testDiffColNum(self):
        diff_list_line = diff_lines(DIFF_LINES_MISS_VAL, DIFF_LINES_BASE_FILE)
        self.assertEqual(len(diff_list_line), 2)

    def testMissLine(self):
        diff_line_list = diff_lines(DIFF_LINES_BASE_FILE, MISS_LINES_MISS_LINE)
        self.assertEqual(len(diff_line_list), 1)
        self.assertTrue("- 540010,1.04337066817119" in diff_line_list[0])

    def testDiffOrd(self):
        diff_line_list = diff_lines(IMPROP_SEC, IMPROP_SEC_ALT, delimiter=" ")
        self.assertEqual(13, len(diff_line_list))

    def testDiffOneNan(self):
        diff_line_list = diff_lines(DIFF_LINES_BASE_FILE, DIFF_LINES_ONE_NAN)
        self.assertEqual(2, len(diff_line_list))

    def testDiffBothNanPrecDiff(self):
        # make there also be a precision difference so the entry-by-entry comparison will be made
        diff_line_list = diff_lines(DIFF_LINES_ONE_NAN_PREC_DIFF, DIFF_LINES_ONE_NAN)
        self.assertFalse(diff_line_list)

    def testStrDiff(self):
        diff = diff_lines(DIFF_LINES_BASE_FILE, DIFF_LINES_STR_DIFF)
        good_diff = ['- 540000,1.0261450032524644,2.23941837370778,1.2132733704553156,1.491028573368929',
                     '+ 540000,1.0261450032524644,2.23941837370778,1.2132733704553156,ghost']
        self.assertEqual(diff, good_diff)

    def testDiffOrder(self):
        # warning("Files differ in trailing white space.")
        test_input = (DIFF_LINES_BASE_FILE, DIFF_LINES_DIFF_ORDER)
        diff_line_list = diff_lines(*test_input)
        self.assertEqual(len(diff_line_list), 2)
        with capture_stderr(diff_lines, *test_input) as output:
            self.assertTrue("order" in output)

    def testDiffWhiteSpace(self):
        temp_file1 = os.path.join(SUB_DATA_DIR, 'temp1.txt')
        temp_file2 = os.path.join(SUB_DATA_DIR, 'temp2.txt')
        test_input = (temp_file1, temp_file2)
        try:
            for fname in test_input:
                with open(fname, "w") as f:
                    f.write("hello")
                    if fname == temp_file1:
                        f.write("\n")
                    else:
                        f.write("     \n")
                    f.write("hello\n")
            diff_line_list = diff_lines(*test_input)
            self.assertEqual(len(diff_line_list), 0)
            with capture_stderr(diff_lines, *test_input) as output:
                self.assertTrue("white space" in output)
        finally:
            for fname in test_input:
                silent_remove(fname, disable=DISABLE_REMOVE)
            pass

    def testDiffPrecisionWhiteSpaceParens(self):
        temp_file1 = os.path.join(SUB_DATA_DIR, 'temp1.txt')
        temp_file2 = os.path.join(SUB_DATA_DIR, 'temp2.txt')
        test_input = (temp_file1, temp_file2)
        try:
            with open(temp_file1, "w") as f:
                f.write("hello\n")
                f.write("(3.121,6.78),(3.121,6.78)\n")
                f.write("\n")
            with open(temp_file2, "w") as f:
                f.write("hello        \n")
                f.write("(3.12099999999999,6.78000000000001),(3.12099999999999,6.78000000000001)\n")
                f.write("\n")
            diff_line_list = diff_lines(*test_input)
            self.assertEqual(len(diff_line_list), 0)
            with capture_stderr(diff_lines, *test_input) as output:
                self.assertTrue("white space" in output)
                self.assertTrue("floating point precision" in output)
        finally:
            for fname in test_input:
                silent_remove(fname, disable=DISABLE_REMOVE)
            pass

    def testDiffPrecisionNan(self):
        temp_file1 = os.path.join(SUB_DATA_DIR, 'temp1.txt')
        temp_file2 = os.path.join(SUB_DATA_DIR, 'temp2.txt')
        test_input = (temp_file1, temp_file2)
        try:
            with open(temp_file1, "w") as f:
                f.write("nan,3.121\n")
                f.write("\n")
            with open(temp_file2, "w") as f:
                f.write("nan,3.12099999999999\n")
                f.write("\n")
            diff_line_list = diff_lines(*test_input)
            self.assertEqual(len(diff_line_list), 0)
            with capture_stderr(diff_lines, *test_input) as output:
                self.assertTrue("floating point precision" in output)
        finally:
            for fname in test_input:
                silent_remove(fname, disable=DISABLE_REMOVE)
            pass

    def testDiffPrecisionStr(self):
        temp_file1 = os.path.join(SUB_DATA_DIR, 'temp1.txt')
        temp_file2 = os.path.join(SUB_DATA_DIR, 'temp2.txt')
        test_input = (temp_file1, temp_file2)
        try:
            with open(temp_file1, "w") as f:
                f.write("fname,3.121\n")
                f.write("\n")
            with open(temp_file2, "w") as f:
                f.write("fname,3.12099999999999\n")
                f.write("\n")
            diff_line_list = diff_lines(*test_input)
            self.assertEqual(len(diff_line_list), 0)
            with capture_stderr(diff_lines, *test_input) as output:
                self.assertTrue("floating point precision" in output)
        finally:
            for fname in test_input:
                silent_remove(fname, disable=DISABLE_REMOVE)
            pass

    def testDiffStr(self):
        temp_file1 = os.path.join(SUB_DATA_DIR, 'temp1.txt')
        temp_file2 = os.path.join(SUB_DATA_DIR, 'temp2.txt')
        test_input = (temp_file1, temp_file2)
        try:
            with open(temp_file1, "w") as f:
                f.write("fname,3000000000\n")
                f.write("\n")
            with open(temp_file2, "w") as f:
                f.write("f_name,3000000000\n")
                f.write("\n")
            diff_line_list = diff_lines(*test_input)
            self.assertEqual(len(diff_line_list), 2)
        finally:
            for fname in test_input:
                silent_remove(fname, disable=DISABLE_REMOVE)
            pass


class TestQuoteDeQuote(unittest.TestCase):
    def testQuoting(self):
        self.assertTrue(quote((0, 1)) == '"(0, 1)"')

    def testNoQuotingNeeded(self):
        self.assertTrue(quote('"(0, 1)"') == '"(0, 1)"')

    def testDequote(self):
        self.assertTrue(dequote('"(0, 1)"') == '(0, 1)')

    def testNoDequoteNeeded(self):
        self.assertTrue(dequote("(0, 1)") == '(0, 1)')

    def testDequoteUnmatched(self):
        self.assertTrue(dequote('"' + '(0, 1)') == '"(0, 1)')

    def testSingleQuote(self):
        self.assertTrue(single_quote("(0, 1)") == "'(0, 1)'")

    def testSingleQuoteAlreadyDone(self):
        self.assertTrue(single_quote("'(0, 1)'") == "'(0, 1)'")

    def testSingleQuoteFromDouble(self):
        self.assertTrue(single_quote('"(0, 1)"') == "'(0, 1)'")


class TestConversions(unittest.TestCase):
    def testNotBool(self):
        try:
            str_to_bool("hello there neighbor")
        except ValueError as e:
            self.assertTrue("Cannot covert" in e.args[0])

    def testIntList(self):
        int_str = '2,3,4'
        int_list = [2, 3, 4]
        self.assertEqual(int_list, conv_raw_val(int_str, []))

    def testNotIntMissFlag(self):
        non_int_str = 'a,b,c'
        try:
            conv_raw_val(non_int_str, [])
        except ValueError as e:
            self.assertTrue("invalid literal for int()" in e.args[0])

    def testNotIntList(self):
        non_int_str = 'a,b,c'
        non_int_list = ['a', 'b', 'c']
        self.assertEqual(non_int_list, conv_raw_val(non_int_str, [], int_list=False))


class TestVectorPBCMath(unittest.TestCase):
    # TODO: consider moving to md_wrangler common
    def testSubtractInSameImage(self):
        self.assertTrue(np.allclose(pbc_calc_vector(VEC_1, VEC_2, PBC_BOX), VEC_21))
        self.assertTrue(np.allclose(pbc_calc_vector(VEC_3, VEC_2, PBC_BOX), VEC_23))
        self.assertFalse(np.allclose(pbc_calc_vector(VEC_3, VEC_2, PBC_BOX), VEC_21))
        self.assertTrue(np.allclose(pbc_calc_vector(A_VEC, B_VEC, PBC_BOX), GOOD_A_MINUS_B))

    def testSubtractInDiffImages(self):
        self.assertTrue(np.allclose(pbc_calc_vector(A_VEC, C_VEC, PBC_BOX), GOOD_A_MINUS_C))

    def testAvgInSameImage(self):
        self.assertTrue(np.allclose(pbc_vector_avg(A_VEC, B_VEC, PBC_BOX), GOOD_A_B_AVG))

    def testAvgInDiffImages(self):
        self.assertTrue(np.allclose(pbc_vector_avg(A_VEC, C_VEC, PBC_BOX), GOOD_A_C_AVG))

    def testUnitVector(self):
        test_unit_vec = unit_vector(VEC_3)
        self.assertTrue(np.allclose(test_unit_vec, UNIT_VEC_3))
        self.assertFalse(np.allclose(test_unit_vec, VEC_1))

    def testAngle(self):
        self.assertAlmostEqual(vec_angle(VEC_21, VEC_23), ANGLE_123)

    def testDihedral(self):
        self.assertAlmostEqual(vec_dihedral(VEC_21, VEC_23, VEC_34), DIH_1234)

    def testDist(self):
        good_dist = 1.2589809371074687
        dist_12 = calc_dist(VEC_1, VEC_2)
        self.assertAlmostEqual(dist_12, good_dist)
        dist_21 = calc_dist(VEC_2, VEC_1)
        self.assertAlmostEqual(dist_21, good_dist)


class TestLongestCommonSubstring(unittest.TestCase):
    def testSameLength(self):
        s1 = "small fur"
        s2 = "Small Fur"
        result = longest_common_substring(s1, s2)
        self.assertTrue(result == "mall ")
        print(result)

    def testDiffLength(self):
        s1 = "small fur"
        s2 = "very small fur"
        result = longest_common_substring(s1, s2)
        self.assertTrue(result == "small fur")

    def testLongerFirst(self):
        s1 = "1 small fur"
        s2 = "very small fur!"
        result = longest_common_substring(s2, s1)
        self.assertTrue(result == " small fur")


class TestReadConfig(unittest.TestCase):
    # Filling it what not covered by fill_tpl
    def testExtraKey(self):
        config = ConfigParser()
        config.read(ONE_KEY_INI)
        def_config_vals = {}
        req_keys = {}
        cfg_dict = process_cfg(dict(config.items(MAIN_SEC)), def_config_vals, req_keys, store_extra_keys=True)
        good_cfg_dict = {'ghost': 'ghost'}
        self.assertTrue(cfg_dict == good_cfg_dict)

    def testMissingReqKey(self):
        try:
            message = ""
            config = ConfigParser()
            config.read(ONE_KEY_INI)
            def_config_vals = {'ghost': 'most'}
            req_keys = {'Learning': bool}
            process_cfg(dict(config.items(MAIN_SEC)), def_config_vals, req_keys)
        except KeyError as e:
            message = e.args[0]
        self.assertTrue("Missing config val for key" in message)

    def testWrongTypeReqKey(self):
        try:
            message = ""
            config = ConfigParser()
            config.read(ONE_KEY_INI)
            def_config_vals = {}
            req_keys = {'ghost': int}
            process_cfg(dict(config.items(MAIN_SEC)), def_config_vals, req_keys)
        except InvalidDataError as e:
            message = e.args[0]
        self.assertTrue("Problem with config vals on key 'ghost'" in message)


class TestChemistry(unittest.TestCase):
    def testStoichCalc1(self):
        stoich1 = 'C6H6'
        good_stoich_dict = {'C': 6, 'H': 6}
        stoich_dict = parse_stoich(stoich1)
        self.assertEqual(stoich_dict, good_stoich_dict)

    def testStoichAddToDict(self):
        ini_stoich_dict = {'C': 6, 'H': 6}
        add_stoich = 'O'
        new_stoich_dict = parse_stoich(add_stoich, add_to_dict=ini_stoich_dict)
        good_stoich_dict = {'O': 1, 'C': 6, 'H': 6}
        self.assertEqual(new_stoich_dict, good_stoich_dict)

    def testStoichAddToDict2(self):
        ini_stoich_dict = {'C': 6, 'H': 6}
        add_stoich = 'CH2'
        new_stoich_dict = parse_stoich(add_stoich, add_to_dict=ini_stoich_dict)
        good_stoich_dict = {'C': 7, 'H': 8}
        self.assertEqual(new_stoich_dict, good_stoich_dict)


class TestAssignColor(unittest.TestCase):
    def testEqualIndex(self):
        for idx in range(0, NUM_COLORS):
            self.assertEqual(assign_color(idx), COLOR_SEQUENCE[idx])

    def testFirstWrap(self):
        for idx in range(0, NUM_COLORS):
            self.assertEqual(assign_color(idx + NUM_COLORS), COLOR_SEQUENCE[idx])

    def testAnotherWrap(self):
        for idx in range(0, NUM_COLORS):
            self.assertEqual(assign_color(idx + NUM_COLORS * 2), COLOR_SEQUENCE[idx])


class TestPlotting(unittest.TestCase):
    """
    make_fig options not yet tested:
        x_fill=None, y_fill=None, x2_fill=None, y2_fill=None,
        fill1_label=None, fill2_label=None,
        fill_color_1="green", fill_color_2="blue",
        fig_width=DEF_FIG_WIDTH, fig_height=DEF_FIG_HEIGHT, axis_font_size=DEF_AXIS_SIZE,
        tick_font_size=DEF_TICK_SIZE, print_msg=True
    Options (with values, not None) that will increase coverage:
        y5_array
        fill1_label, fill2_label
        x_fill, x2_fill
    """

    def testSimplePlot(self):
        # Smoke test only, that there are no errors using these options
        silent_remove(TEST_PLOT_FNAME)
        x_values = [0.02, 1.27, 2.28, 3.27, 4.3, 5.32, 6.35, 7.4, 8.08, 9.71, 10.69]
        y_values = [75901, 616236, 304071, 880000, 864863, 299723, 297125, 289680, 684642, 256205, 1320600]
        x_label = "Time (min)"
        y_label = "Intensities (unscaled)"
        try:
            make_fig(TEST_PLOT_FNAME, x_values, y_values, x_label=x_label, y_label=y_label, loc=0)
            self.assertTrue(os.path.isfile(TEST_PLOT_FNAME))
        finally:
            silent_remove(TEST_PLOT_FNAME, disable=DISABLE_REMOVE)

    def testPlotMultYVals(self):
        # Smoke test only, that there are no errors using these options
        silent_remove(TEST_PLOT_FNAME)
        title = 'Free Energy Diagram'
        x_values = np.array([0,  1, 3, 4, 6, 7, 9, 10, 12, 13])
        y_values = [np.array([0., 0., 41.16, 41.16, -4.54, -4.54, 39.72, 39.72, -5.61, -5.62]),
                    np.array([0., 0., 42.76, 42.75, -3.42, -3.42, 40.39, 40.39, -5.80, -5.80]),
                    np.array([0., 0., 41.13, 41.14, -1.05, -1.05, 40.16, 40.16, -2.26, -2.26]),
                    np.array([0., 0., 37.44, 37.44, -4.33, -4.33, 33.87, 33.87, -11.38, -11.38]),
                    None]
        x_label = 'reaction coordinate'
        y_label = '\u0394G at 460 K (kcal/mol)'
        y_labels = ['TPA', 'IPA', 'PDC1', 'PDC2', None]
        try:
            make_fig(TEST_PLOT_FNAME, x_values, y_values[0], x_label=x_label, y_label=y_label,
                     y1_label=y_labels[0], y2_label=y_labels[1], y3_label=y_labels[2], y4_label=y_labels[3],
                     y5_label=y_labels[4], y2_array=y_values[1], y3_array=y_values[2], y4_array=y_values[3],
                     y5_array=y_values[4], ls2='-', ls3='-', ls4='-', ls5='-',
                     x_limb=14, y_lima=-15, y_limb=45, loc=0, hide_x=True, title=title)
            self.assertTrue(os.path.isfile(TEST_PLOT_FNAME))
        finally:
            silent_remove(TEST_PLOT_FNAME, disable=DISABLE_REMOVE)
            pass

    def testPlotIncreaseCoverage(self):
        # Smoke test only, that there are no errors using options: y5_array, y_limb without y_lima
        silent_remove(TEST_PLOT_FNAME)
        title = 'Free Energy Diagram'
        x_values = np.array([0,  1, 3, 4, 6, 7, 9, 10, 12, 13])
        y_values = [np.array([0., 0., 41.16, 41.16, -4.54, -4.54, 39.72, 39.72, -5.61, -5.62]),
                    np.array([0., 0., 42.76, 42.75, -3.42, -3.42, 40.39, 40.39, -5.80, -5.80]),
                    np.array([0., 0., 41.13, 41.14, -1.05, -1.05, 40.16, 40.16, -2.26, -2.26]),
                    np.array([0., 0., 37.44, 37.44, -4.33, -4.33, 33.87, 33.87, -11.38, -11.38]),
                    np.array([0., 0., 37.44, 37.44, -4.33, -4.33, 13.87, 13.87, -10.38, -10.38])]
        x_label = 'reaction coordinate'
        y_label = '\u0394G at 460 K (kcal/mol)'
        y_labels = ['TPA', 'IPA', 'PDC1', 'PDC2', None]
        try:
            make_fig(TEST_PLOT_FNAME, x_values, y_values[0], x_label=x_label, y_label=y_label,
                     y1_label=y_labels[0], y2_label=y_labels[1], y3_label=y_labels[2], y4_label=y_labels[3],
                     y5_label=y_labels[4], y2_array=y_values[1], y3_array=y_values[2], y4_array=y_values[3],
                     y5_array=y_values[4], ls2='-', ls3='-', ls4='-', ls5='-',
                     x_limb=14, y_limb=45, loc=0, hide_x=True, title=title)
            self.assertTrue(os.path.isfile(TEST_PLOT_FNAME))
        finally:
            silent_remove(TEST_PLOT_FNAME, disable=DISABLE_REMOVE)
            pass


class TestConvStrToFunc(unittest.TestCase):
    def testToSTR(self):
        input_str = 'str'
        test_val = 123.1
        good_out = str(test_val)
        # function to test on next line
        out_func = conv_str_to_func(input_str)
        out_val = out_func(test_val)
        self.assertEqual(out_val, good_out)
        # next to make sure the test was valid...
        self.assertNotEqual(test_val, out_val)

    def testToInt(self):
        input_str = 'int'
        test_val = '124'
        good_out = int(test_val)
        # function to test on next line
        out_func = conv_str_to_func(input_str)
        out_val = out_func(test_val)
        self.assertEqual(out_val, good_out)
        self.assertNotEqual(test_val, out_val)

    def testToFloat(self):
        input_str = 'float'
        test_val = '123.1'
        good_out = float(test_val)
        # function to test on next line
        out_func = conv_str_to_func(input_str)
        out_val = out_func(test_val)
        self.assertAlmostEqual(out_val, good_out)
        self.assertNotEqual(test_val, out_val)

    def testToBoolean(self):
        input_str = 'bool'
        test_val = 'False'
        good_out = bool(test_val)
        # function to test on next line
        out_func = conv_str_to_func(input_str)
        out_val = out_func(test_val)
        self.assertEqual(out_val, good_out)
        self.assertNotEqual(test_val, out_val)

    def testNoneStr(self):
        input_str = 'None'
        out_func = conv_str_to_func(input_str)
        self.assertIsNone(out_func)

    def testNone(self):
        input_str = None
        out_func = conv_str_to_func(input_str)
        self.assertIsNone(out_func)

    def testNonValidInput(self):
        input_str = 'ghost'
        try:
            # function to test on next line
            conv_str_to_func(input_str)
            self.assertFalse("I should not get here...")
        except InvalidDataError as e:
            error_str = e.args[0]
            self.assertTrue(input_str in error_str)
