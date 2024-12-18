#!/usr/bin/env python
# -*- coding: utf-8 -*-


from sequenza.misc import get_modules, DefaultHelpParser, SubcommandHelpFormatter
from sequenza import __version__
from sequenza.misc import SeqzLogger
import sequenza.programs


def main():
    """
    Execute the function with args
    """
    parser = DefaultHelpParser(
        prog="sequenza-utils",
        formatter_class=lambda prog: SubcommandHelpFormatter(
            prog, max_help_position=20, width=75
        ),
        description=(
            "Sequenza Utils is a collection of tools primarily "
            "design to convert bam, pileup and vcf files to seqz "
            "files, the format used in the sequenza R package"
        ),
        add_help=True,
        epilog="This is version %s - %s - %s"
        % (__version__.VERSION, __version__.AUTHOR, __version__.DATE),
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        help="Show all logging information",
        action="store_true",
    )
    subparsers = parser.add_subparsers(dest="module")

    modules = get_modules(sequenza.programs, subparsers, {})
    try:
        args, extra = parser.parse_known_args()
        if args.verbose is True:
            log_level = 1
        else:
            log_level = 30
        if args.module in modules.keys():
            log = SeqzLogger(level=log_level)
            modules[args.module](subparsers, args.module, extra, log)
        else:
            if args.module is None:
                return parser.print_help()
            else:
                return parser.parse_args(args)
    except IndexError:
        return parser.print_help()


if __name__ == "__main__":
    main()
