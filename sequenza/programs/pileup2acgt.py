from sequenza.misc import xopen

try:
    import __pypy__
    from sequenza.pileup import acgt
except ImportError:
    try:
        from sequenza.c_pileup import acgt
    except ImportError:
        from sequenza.pileup import acgt


def add_parser(subparsers, module_name):
    return subparsers.add_parser(
        module_name,
        add_help=False,
        help=(
            "Parse the format from the samtools "
            "mpileup command, and report the "
            "occurrence of the 4 nucleotides "
            "in each position."
        ),
    )


def pileup2actg_args(parser, argv):
    parser_io = parser.add_argument_group(
        title="Output", description="Arguments that involve the output destination."
    )
    parser_io.add_argument(
        "-p",
        "--mpileup",
        dest="mpileup",
        required=True,
        help=(
            "Name of the input mpileup (SAMtools) file. If the "
            "filename ends in .gz it will be opened in gzip mode. If "
            "the file name is - it will be read from STDIN."
        ),
    )
    parser_io.add_argument(
        "-o",
        "--output",
        dest="output",
        default="-",
        help=(
            "Name of the output file. To use gzip compression name the "
            "file ending in .gz. Default STDOUT."
        ),
    )
    parser_filters = parser.add_argument_group(
        title="Filtering and Format",
        description="Arguments that apply some filter to process the mpileup.",
    )
    parser_filters.add_argument(
        "-n",
        dest="n",
        default=1,
        type=int,
        help=("The minimum required read depth on a position to test " "for mutation."),
    )
    parser_filters.add_argument(
        "-q",
        "--qlimit",
        dest="qlimit",
        default=20,
        type=int,
        help="Minimum nucleotide quality score filter.",
    )
    parser_filters.add_argument(
        "--no-end",
        dest="noend",
        action="store_true",
        help="Discard the base located at the end of the read",
    )
    parser_filters.add_argument(
        "--no-start",
        dest="nostart",
        action="store_true",
        help="Discard the base located at the start of the read",
    )
    parser_filters.add_argument(
        "-f",
        "--qformat",
        dest="qformat",
        default="sanger",
        help=(
            'Quality format, options are "sanger" or "illumina". '
            "This will add an offset of 33 or 64 respectively to "
            "the qlimit value."
        ),
    )

    return parser.parse_args(argv)


def pileup2acgt(subparsers, module_name, argv, log):
    log.log.info("Start %s" % module_name)
    args = pileup2actg_args(add_parser(subparsers, module_name), argv)
    header = ["chr", "n_base", "ref_base", "read.depth", "A", "C", "G", "T", "strand"]
    qlimit = args.qlimit
    if args.qformat == "sanger":
        qlimit = qlimit + 33
    elif args.qformat == "illumina":
        qlimit = qlimit + 64
    else:
        raise Exception(
            ('Supported quality format are "illumina" ' 'and "sanger"(default).')
        )
    acgt_format = "%s\t%s\t%s\t%i\t%i\t%i\t%i\t%i\t%i:%i:%i:%i\n"
    with xopen(args.mpileup, "rt") as mpileup, xopen(args.output, "wt") as out:
        out.write("%s\n" % "\t".join(map(str, header)))
        for line in mpileup:
            try:
                (
                    chromosome,
                    position,
                    reference,
                    depth,
                    pileup,
                    quality,
                ) = line.strip().split("\t")
                depth = int(depth)
                reference = reference.upper()
                if depth >= args.n and args.n > 0 and reference != "N":
                    acgt_res = acgt(
                        pileup,
                        quality,
                        depth,
                        reference,
                        qlimit=qlimit,
                        noend=args.noend,
                        nostart=args.nostart,
                    )
                    out.write(
                        acgt_format
                        % (
                            chromosome,
                            position,
                            reference,
                            depth,
                            acgt_res["A"],
                            acgt_res["C"],
                            acgt_res["G"],
                            acgt_res["T"],
                            acgt_res["Z"][0],
                            acgt_res["Z"][1],
                            acgt_res["Z"][2],
                            acgt_res["Z"][3],
                        )
                    )
            except ValueError:
                pass
