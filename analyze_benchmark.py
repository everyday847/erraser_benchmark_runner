"""
ERRASER benchmark analysis, v. 0.1
The idea is, you hand it a folder with a PDB and an MTZ. It figures out your state,
it cleans up your messes, it does analysis.
"""

import argparse
from logging import info, warning, error, basicConfig, INFO
import os
import glob
import functools # for partial
import shutil # for rmtree
from typing import Callable, List

from molprobity import Result, Clash

def analyze(d: str):
    """
    Create Result objects from the molprobity files in the provided directory.
    Make all before/after comparisons.
    """

    starting_molprobity = Result("{}/{}_start.molprobity".format(d, d))
    # will have to supply nstruct or figure out the range here
    end_molprobity = { ii: Result("{}/output/runs/{}/{}_end.molprobity".format(d, ii, d)) for ii in range(100) }

    for ii, end_mol in end_molprobity.items():
        for c1, c2 in zip(starting_molprobity.clashes, end_mol.clashes):
            pass
        print("{} -- {} => {}".format(ii, starting_molprobity.clashscore, end_mol.clashscore))

def main(directories: List[str]) -> None:
    for d in directories:
        info("Analyzing ERRASER2 on contents of directory {}.".format(d))
        analyze(d)

def test_args():
    pass

if __name__ == '__main__':

    basicConfig(level=INFO)

    parser = argparse.ArgumentParser(description='Run ERRASER2 on a working directory')
    parser.add_argument('directory', type=str, nargs='+',
                        help='dir(s) to run')
    args = parser.parse_args()
    main(args.directory)
