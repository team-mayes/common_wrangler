# coding=utf-8

"""
"""
import os
import unittest
import numpy as np
import matplotlib.pyplot as plt
from common_wrangler.common import COLORBREWER_BLUE, COLORBREWER_GREEN, COLORBREWER_PURPLE, make_dir

__author__ = 'hmayes'


# MAKE_PLOTS = True
MAKE_PLOTS = False


# Directories #

DATA_DIR = os.path.join(os.path.dirname(__file__), 'test_data')
SUB_DATA_DIR = os.path.join(DATA_DIR, 'scaling')

# Files #


# Common methods #

def plot_scaling(nodes, ns_per_day, x_label, x_ticks, plot_f_loc, location='upper center'):
    ideal_ns_per_day = nodes * ns_per_day[0]
    efficiency = ns_per_day / ideal_ns_per_day
    color1 = COLORBREWER_BLUE
    color2 = COLORBREWER_GREEN
    color3 = COLORBREWER_PURPLE

    fig, ax1 = plt.subplots(figsize=(5.0, 2.8))
    ax1.set_xlabel(x_label)
    ax1.set_ylabel('ns/day', color=color1)
    ax1.plot(nodes, ns_per_day, color=color1, marker='o', label="real scaling")
    ax1.plot(nodes, ideal_ns_per_day, color=color2, label="ideal scaling")
    ax1.tick_params(axis='y', labelcolor=color1)
    ax1.set_ylim([0., ideal_ns_per_day[-1]])

    ax2 = ax1.twinx()
    ax2.set_ylabel('efficiency (% real vs. ideal)', color=color3)
    ax2.plot(nodes, efficiency, color=color3, marker='o', label="efficiency")
    ax2.tick_params(axis='y', labelcolor=color3)
    ax2.set_ylim([0., 1.])

    plt.xticks(x_ticks)
    plt.xlim([nodes[0], nodes[-1]])

    if location:
        fig.legend(loc=location)
    fig.tight_layout()
    plt.savefig(plot_f_loc, format='png', transparent=True, bbox_inches='tight')
    plt.close()


# Tests #

class TestPlotData(unittest.TestCase):
    def testCel7A(self):
        if MAKE_PLOTS:
            make_dir(SUB_DATA_DIR)
            # Amber16, ~70,000 atoms
            nodes = np.asarray([1, 2, 4, 8, 16, 32, 64, 128, 256, 512])
            ns_per_day = np.asarray([0.14, 0.28, 0.55, 1.07, 2.11, 3.97, 6.38, 9.49, 14.20, 13.16])
            x_label = 'Number of Comet Cores'
            x_ticks = [1, 36, 64, 128, 256, 512]
            plot_f_loc = os.path.join(SUB_DATA_DIR, "amber_md_cel7a.png")
            plot_scaling(nodes, ns_per_day, x_label, x_ticks, plot_f_loc)

    def testCESA(self):
        if MAKE_PLOTS:
            make_dir(SUB_DATA_DIR)
            # GROMACS, ~128,000 atoms
            nodes = np.asarray([1, 2, 4, 8])
            ns_per_day = np.asarray([10.6, 20.7, 40.1, 55.7])
            x_label = 'Number of Bridges Nodes'
            x_ticks = [1, 2, 4, 8]
            plot_f_loc = os.path.join(SUB_DATA_DIR, "gromacs_md_cesa.png")
            plot_scaling(nodes, ns_per_day, x_label, x_ticks, plot_f_loc, location="")

    def testCelluloseFibril(self):
        if MAKE_PLOTS:
            make_dir(SUB_DATA_DIR)
            # NAMD, ~1.3 million atoms
            nodes = np.asarray([1, 2, 4, 8, 16])
            ns_per_day = np.asarray([0.34, 0.62, 1.22, 2.38, 4.71])
            x_label = 'Number of Bridges Nodes'
            x_ticks = nodes
            plot_f_loc = os.path.join(SUB_DATA_DIR, "gromacs_md_fibril.png")
            plot_scaling(nodes, ns_per_day, x_label, x_ticks, plot_f_loc, location="")

    def testTmAfc(self):
        if MAKE_PLOTS:
            make_dir(SUB_DATA_DIR)
            cores = np.array([1, 2, 4, 8, 16, 24, 48])
            ns_per_day = np.array([0.5, 1, 1.9, 3.65, 6.57, 9.36, 9.4])
            x_label = 'Number of Comet Cores'
            x_ticks = cores
            plot_f_loc = os.path.join(SUB_DATA_DIR, "amber_tmafc.png")
            plot_scaling(cores, ns_per_day, x_label, x_ticks, plot_f_loc, location="")

    def testDFT(self):
        if MAKE_PLOTS:
            make_dir(SUB_DATA_DIR)
            cores = np.array([7, 14, 28])
            min_per_sp = np.array([39.25166667, 21.50833333, 11.55833333])
            x_label = 'Number of Bridges Cores'
            x_ticks = cores
            y1_label = 'Minutes per single-point calculation'
            plot_f_loc = os.path.join(SUB_DATA_DIR, "bridges_dft.png")

            ideal_cores = np.array([1, 7, 14, 28])
            ideal_min_per_sp = np.array([274.7616667, 39.25166667, 19.62583333, 9.812916667])
            efficiency = np.array([1., 1., .912, .848990627])
            color1 = COLORBREWER_BLUE
            color2 = COLORBREWER_GREEN
            color3 = COLORBREWER_PURPLE

            fig, ax1 = plt.subplots(figsize=(5.0, 2.8))
            ax1.set_xlabel(x_label)
            ax1.set_ylabel(y1_label, color=color1)
            ax1.plot(cores, min_per_sp, color=color1, marker='o', label="real scaling")
            ax1.plot(ideal_cores, ideal_min_per_sp, color=color2, label="ideal scaling")
            ax1.tick_params(axis='y', labelcolor=color1)
            ax1.set_ylim([0., 40])

            ax2 = ax1.twinx()
            ax2.set_ylabel('efficiency (% real vs. ideal)', color=color3)
            ax2.plot(ideal_cores, efficiency, color=color3, marker='o', label="efficiency")
            ax2.tick_params(axis='y', labelcolor=color3)
            ax2.set_ylim([0., 1.])

            plt.xticks(x_ticks)
            plt.xlim([1, 28])

            fig.tight_layout()
            plt.savefig(plot_f_loc, format='png', transparent=True, bbox_inches='tight')
            plt.close()
