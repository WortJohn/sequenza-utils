from sequenza.misc import xopen
from sequenza.seqz import binned_seqz
from sequenza.samtools import tabix_seqz


def add_parser(subparsers, module_name):
    return subparsers.add_parser(
        module_name,
        add_help=False,
        help=(
            "Perform the binning of the seqz file to reduce file size "
            "and memory requirement for the analysis."
        ),
    )


def binning_args(parser, argv):
    parser.add_argument("-s", "--seqz", dest="seqz", required=True, help="A seqz file.")
    parser.add_argument(
        "-w",
        "--window",
        dest="window",
        type=int,
        default=50,
        help=("Window size used for binning the original " "seqz file. Default is 50."),
    )
    parser.add_argument(
        "-o", dest="out", type=str, default="-", help='Output file "-" for STDOUT'
    )

    parser.add_argument(
        "-T",
        "--tabix",
        dest="tabix",
        type=str,
        default="tabix",
        help=("Path of the tabix binary. " 'Default "tabix"'),
    )
    return parser.parse_args(argv)


def seqz_binning(subparsers, module_name, argv, log):
    log.log.info("Start %s" % module_name)
    args = binning_args(add_parser(subparsers, module_name), argv)
    with xopen(args.seqz, "rt") as seqz, xopen(args.out, "wt", bgzip=True) as out:
        bins = binned_seqz(seqz, args.window)
        out.write(next(bins))
        for line in bins:
            out.write(line)
    if args.out.endswith(".gz"):
        tabix_seqz(args.out, args.tabix)
