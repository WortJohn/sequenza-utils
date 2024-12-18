from sequenza.misc import xopen
from sequenza.vcf import vfc2seqz, vcf_parse
from sequenza.samtools import tabix_seqz
from sequenza.wig import Wiggle


def add_parser(subparsers, module_name):
    return subparsers.add_parser(
        module_name,
        add_help=False,
        help=(
            "Parse VCFs and other variant and "
            "coverage formats to produce seqz "
            "files"
        ),
    )


def snp2seqz_args(parser, argv):
    parser_io = parser.add_argument_group(
        title="Output", description="Output arguments"
    )
    parser_in = parser.add_argument_group(title="Input", description="Input files")
    parser_vcf = parser.add_argument_group(
        title="VCF", description="Parsing option for the VCF file"
    )
    parser_genotype = parser.add_argument_group(
        title="Genotype", description="Genotype filtering options"
    )
    parser_programs = parser.add_argument_group(
        title="Programs", description="Option to call and control externa programs"
    )
    parser_filter = parser.add_argument_group(
        title="Filters", description="Filter output file by various parameters"
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
        "-v", "--vcf", dest="vcf", required=True, help=("VCF input file")
    )
    parser_in.add_argument(
        "-gc", dest="gc", required=True, help="The GC-content wiggle file"
    )
    parser_vcf.add_argument(
        "--vcf-depth",
        dest="vcf_depth_tag",
        type=str,
        default="DP:DP",
        help=(
            "Column separated VCF tags in the format column "
            "for the read depth for the normal and for the tumor. "
            'Default "DP:DP"'
        ),
    )
    parser_vcf.add_argument(
        "--vcf-samples",
        dest="vcf_samples_order",
        type=str,
        choices=["n/t", "t/n"],
        default="n/t",
        help=(
            "Order of the normal and tumor sample in the VCF column, "
            'choices are "n/t" or "t/n". Default "n/t"'
        ),
    )
    parser_vcf.add_argument(
        "--vcf-alleles",
        dest="vcf_alleles_tag",
        type=str,
        default="AD:AD",
        help=(
            "Column separated VCF tags in the format column "
            "for the alleles depth for the normal and for the tumor. "
            'Default "AD:AD"'
        ),
    )
    parser_vcf.add_argument(
        "--preset",
        dest="vcf_preset",
        type=str,
        choices=["caveman", "mutect", "mpileup", "strelka2_som", "dp4"],
        help="Preset set of options to parse VCF from popular variant callers",
        default=None,
    )
    parser_genotype.add_argument(
        "--hom",
        dest="hom",
        type=float,
        default=0.9,
        help="Threshold to select homozygous positions. Default 0.9",
    )
    parser_genotype.add_argument(
        "--het",
        dest="het",
        type=float,
        default=0.25,
        help="Threshold to select heterozygous positions. Default 0.25.",
    )
    parser_genotype.add_argument(
        "--het_f",
        dest="het_f",
        type=float,
        default=-0.2,
        help=(
            "Threshold of frequency in the forward strand to trust "
            "heterozygous calls. Default -0.2 (Disabled, "
            "effective with values >= 0)."
        ),
    )
    parser_genotype.add_argument(
        "--exclude_hom",
        dest="het_only",
        action="store_true",
        help=("Exclude mutations and keep only heterozygous sites."),
    )
    parser_filter.add_argument(
        "-N",
        dest="n",
        type=int,
        default=20,
        help=(
            "Threshold to filter positions by the sum of read "
            "depth of the two samples. Default 20."
        ),
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


def snp2seqz(subparsers, module_name, argv, log):
    log.log.info("Start %s" % module_name)
    subp = add_parser(subparsers, module_name)
    args = snp2seqz_args(subp, argv)
    with xopen(args.output, "wt", bgzip=True) as out, xopen(
        args.vcf, "rt"
    ) as vcf_input, xopen(args.gc, "rt") as gc:
        vcf_file = vcf_parse(
            vcf_input,
            args.vcf_samples_order,
            "FORMAT",
            args.vcf_depth_tag.split(":"),
            args.vcf_alleles_tag.split(":"),
            args.vcf_preset,
        )
        gc_wig = Wiggle(gc)
        seqz_vcf = vfc2seqz(
            vcf_file, gc_wig, args.hom, args.het, args.het_f, args.het_only
        )
        out.write("%s\n" % "\t".join(next(seqz_vcf)))
        for vcf_line in seqz_vcf:
            if vcf_line[3] + vcf_line[4] >= args.n:
                out.write("%s\n" % "\t".join(map(str, vcf_line)))
    if args.output.endswith(".gz"):
        tabix_seqz(args.output, args.tabix)
