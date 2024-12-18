from sequenza.misc import xopen, split_coordinates
from sequenza.samtools import tabix_seqz
from sequenza.izip import zip_fast


def add_parser(subparsers, module_name):
    return subparsers.add_parser(
        module_name, add_help=False, help=("Merge two seqz files")
    )


def seqz_merge_args(parser, argv):
    parser_io = parser.add_argument_group(
        title="Output", description="Output arguments"
    )
    parser_in = parser.add_argument_group(title="Input", description="Input files")
    parser_programs = parser.add_argument_group(
        title="Programs", description="Option to call and control externa programs"
    )
    parser_io.add_argument(
        "-o",
        "--output",
        dest="output",
        default="-",
        help=(
            "Output file. For gzip compressed output name the "
            "file ending in .gz. Default STDOUT"
        ),
    )
    parser_in.add_argument(
        "-1",
        "--seqz1",
        dest="s1",
        required=True,
        help=(
            "First input file. If both input files "
            "contain the same line, the information "
            "in the first file will be used"
        ),
    )

    parser_in.add_argument(
        "-2", "--seqz2", dest="s2", required=True, help=("Second input file")
    )
    parser_programs.add_argument(
        "-T",
        "--tabix",
        dest="tabix",
        type=str,
        default="tabix",
        help=('Path of the tabix binary. Default "tabix"'),
    )

    return parser.parse_args(argv)


def seqz_merge(subparsers, module_name, argv, log):
    log.log.info("Start %s" % module_name)
    subp = add_parser(subparsers, module_name)
    args = seqz_merge_args(subp, argv)
    with xopen(args.output, "wt", bgzip=True) as out, xopen(
        args.s1, "rt"
    ) as input_1, xopen(args.s2, "rt") as input_2:
        header = next(input_1)
        header = next(input_2)
        input_1s = split_coordinates(input_1)
        input_2s = split_coordinates(input_2)
        seqz_merged = zip_fast(input_1s, input_2s)
        out.write(header)
        for line in seqz_merged:
            coord = (line[0][0], line[0][1])
            data = (line[1][0],)
            out.write("%s\n" % "\t".join(map(str, coord + data)))
    if args.output.endswith(".gz"):
        tabix_seqz(args.output, args.tabix)
