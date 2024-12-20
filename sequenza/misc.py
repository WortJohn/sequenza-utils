from __future__ import division

import os
import sys
import gzip
from sequenza.bgzf import BgzfWriter
from sequenza import args as argparse
import logging


def xopen(filename, mode="r", bgzip=False):
    """
    Replacement for the "open" function that can also open
    files that have been compressed with gzip. If the filename ends with .gz,
    the file is opened with gzip.open(). If it doesn't, the regular open()
    is used. If the filename is '-', standard output (mode 'w') or input
    (mode 'r') is returned.
    """
    assert isinstance(filename, str)
    if filename == "-":
        return sys.stdin if "r" in mode else sys.stdout
    if filename.endswith(".gz"):
        if bgzip and mode.startswith("w"):
            return BgzfWriter(filename, mode)
        else:
            return gzip.open(filename, mode)
    else:
        return open(filename, mode)


def split_ext(path):
    name, extension = os.path.splitext(path)
    if extension in [".bz2", ".gz", ".txt"]:
        name, extension2 = os.path.splitext(name)
        extension = extension2 + extension
    return name, extension


def split_coordinates(file_in):
    """
    Generator to split a file that have chromosome and position
    as first two columns, keep the remaining part of the line
    as it is.
    The resulting object is a tuple: ((chromosome, position), (data, ))
    """
    if file_in is not None:
        for line in file_in:
            chromosome, position, data = line.strip().split("\t", 2)
            position = int(position)
            yield ((chromosome, position, position + 1), (data,))


def countN(x, n):
    gc = x.count(n)
    return 100 * gc / len(x)


def package_modules(package):
    pathname = package.__path__[0]
    return set(
        [
            ".".join([package.__name__, os.path.splitext(module)[0]])
            for module in os.listdir(pathname)
            if module.endswith(".py") and not module.startswith("__init__")
        ]
    )


def get_modules(parent, subparsers, progs):
    """
    return the list of modules in the program module
    """
    mods = package_modules(parent)
    for mod in sorted(mods):
        try:
            __import__(mod)
            mod_name = mod.split(".")[-1]
            m = getattr(parent, mod_name)
            m.add_parser(subparsers, mod_name)
            progs[mod_name] = getattr(m, mod_name)
        except AttributeError:
            pass
    return progs


class SeqzLogger:
    """
    Custom logger for seuenza utils
    """

    def __init__(self, level):
        self.__name__ = "sequenza-utils"
        self.log = logging.getLogger(self.__name__)
        formatter = logging.Formatter("%(levelname)s : %(asctime)s : %(message)s")
        streamHandler = logging.StreamHandler()
        streamHandler.setFormatter(formatter)
        self.log.setLevel(level)
        self.log.addHandler(streamHandler)


class SubcommandHelpFormatter(argparse.RawDescriptionHelpFormatter):
    def _format_action(self, action):
        parts = super(argparse.RawDescriptionHelpFormatter, self)._format_action(action)
        if action.nargs == argparse.PARSER:
            parts = "\n".join(parts.split("\n")[1:])
        return parts


class DefaultHelpParser(argparse.ArgumentParser):
    def error(self, message):
        import sys

        sys.stderr.write("error: %s\n" % message)
        self.print_help()
        sys.exit(2)
