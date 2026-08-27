#  This source file is part of the Avogadro project
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").
#  https://github.com/ghutchis/avogadro-stk/

"""Supramolecular construction for Avogadro 2, using stk."""

import argparse

from .run import run


def main():
    # Avogadro calls the plugin as:
    #   avogadro-stk <identifier> [--lang <locale>] [--debug]
    # with the user options JSON on stdin.
    parser = argparse.ArgumentParser("avogadro-stk")
    parser.add_argument("feature")
    parser.add_argument("--lang", nargs="?", default="en")
    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()
    run(args)
