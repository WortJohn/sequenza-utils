from sequenza.misc import xopen
from sequenza.bed import SimpleBed, seqz_from_beds
from sequenza.samtools import tabix_seqz
from sequenza.seqz import seqz_header
from sequenza.wig import Wiggle
from sequenza.izip import zip_coordinates


def add_parser(subparsers, module_name):
    return subparsers.add_parser(
        module_name,
        add_help=False,
        help=(
            "Parse BED file from mosdepth (and others) "
            "to produce a seqz file with basic coverage "
            "information"
        ),
    )


def mosdepth2seqz_args(parser, argv):
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
        "-n", "--normal", dest="normal", required=True, help=("Mosdepth BED file for the normal sample")
    )
    parser_in.add_argument(
        "-t", "--tumor", dest="tumor", required=True, help=("Mosdepth BED file for the tumor sample")
    )
    parser_in.add_argument(
        "-gc", dest="gc", required=True, help="The GC-content wiggle file"
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


def mosdepth2seqz(subparsers, module_name, argv, log):
    log.log.info("Start %s" % module_name)
    subp = add_parser(subparsers, module_name)
    args = mosdepth2seqz_args(subp, argv)
    with xopen(args.output, "wt", bgzip=True) as out, xopen(
        args.normal, "rt"
    ) as normal, xopen(args.tumor, "rt") as tumor, xopen(args.gc, "rt") as gc:
        gc_wig = Wiggle(gc)
        bed_normal = SimpleBed(normal)
        bed_tumor = SimpleBed(tumor)
        bed_covs = zip_coordinates(bed_normal, bed_tumor)
        bed_gc = zip_coordinates(bed_covs, gc_wig)
        out.write("%s\n" % "\t".join(seqz_header()))
        for bed_line in seqz_from_beds(bed_gc):
            out.write(bed_line)
    if args.output.endswith(".gz"):
        tabix_seqz(args.output, args.tabix)
