from sequenza.misc import xopen, split_ext, split_coordinates
from sequenza.izip import zip_coordinates
from sequenza.wig import Wiggle
from sequenza.seqz import do_seqz, seqz_header
from sequenza.samtools import bam_mpileup, indexed_pileup, tabix_seqz
from multiprocessing import Pool
from functools import partial
from contextlib import closing


def add_parser(subparsers, module_name):
    return subparsers.add_parser(
        module_name,
        add_help=False,
        help=(
            "Process a paired set of BAM/pileup files (tumor and "
            "matching normal), and GC-content genome-wide "
            "information, to extract the common positions with"
            "A and B alleles frequencies"
        ),
    )


def bam2seqz_args(parser, argv):
    parser_input = parser.add_argument_group(
        title="Input/Output", description="Input and output files."
    )
    parser_genotype = parser.add_argument_group(
        title="Genotype", description="Options regarding the genotype filtering."
    )
    parser_subset = parser.add_argument_group(
        title="Subset indexed files",
        description="Option regarding samtools and bam indexes.",
    )
    parser_qualitysets = parser.add_argument_group(
        title="Quality and Format",
        description="Options that change the quality threshold and format.",
    )
    parser_input.add_argument(
        "-p",
        "--pileup",
        dest="pileup",
        action="store_true",
        help="Use pileups as input files instead of BAMs.",
    )
    parser_input.add_argument(
        "-n",
        "--normal",
        dest="normal",
        required=True,
        help="Name of the BAM/pileup file from the reference/normal sample",
    )
    parser_input.add_argument(
        "-t",
        "--tumor",
        dest="tumor",
        required=True,
        help="Name of the BAM/pileup file from the tumor sample",
    )
    parser_input.add_argument(
        "-gc", dest="gc", required=True, help="The GC-content wiggle file"
    )
    parser_input.add_argument(
        "-F",
        "--fasta",
        dest="fasta",
        default=None,
        help=(
            "The reference FASTA file used to generate the intermediate "
            "pileup. Required when input are BAM"
        ),
    )

    parser_input.add_argument(
        "-o",
        "--output",
        dest="out",
        default="-",
        help=(
            "Name of the output file. To use gzip compression name the "
            "file ending in .gz. Default STDOUT."
        ),
    )
    parser_input.add_argument(
        "-n2",
        "--normal2",
        dest="normal2",
        type=str,
        default=None,
        help=(
            "Optional BAM/pileup file used to compute the "
            "depth.normal and depth-ratio, instead of "
            "using the normal BAM."
        ),
    )
    parser_subset.add_argument(
        "-C",
        "--chromosome",
        dest="chr",
        nargs="+",
        default=[],
        help=(
            "Argument to restrict the input/output to a chromosome or a "
            "chromosome region. Coordinate format is "
            "Name:pos.start-pos.end, eg: chr17:7565097-7590856, for a "
            "particular region; eg: chr17, for the entire chromosome. "
            "Chromosome names can checked in the BAM/pileup files and "
            "are depending on the FASTA reference used for alignment. "
            "Default behavior is to not selecting any chromosome."
        ),
    )
    parser_subset.add_argument(
        "-L", 
        "--target_file",
        dest="target_file",
        default=None,
        help=(
            "overlap (BED) regions in TARGET_FILE"
            ),
        )
    parser_subset.add_argument(
        "--parallel",
        dest="nproc",
        type=int,
        default=1,
        help=(
            "Defines the number of chromosomes to run in parallel. "
            "The output will be divided in multiple files, one for "
            "each chromosome. The file name will be composed by the output "
            "argument (used as prefix) and a chromosome name given by the "
            "chromosome argument list. This imply that both output and "
            "chromosome argument need to be set correctly."
        ),
    )
    parser_subset.add_argument(
        "-S",
        "--samtools",
        dest="samtools",
        type=str,
        default="samtools",
        help=(
            "Path of samtools exec file to access the indexes and compute"
            ' the pileups. Default "samtools"'
        ),
    )
    parser_subset.add_argument(
        "-T",
        "--tabix",
        dest="tabix",
        type=str,
        default="tabix",
        help=("Path of the tabix binary. " 'Default "tabix"'),
    )
    parser_qualitysets.add_argument(
        "-q",
        "--qlimit",
        dest="qlimit",
        default=20,
        type=int,
        help=(
            "Minimum nucleotide quality score for inclusion in the "
            "counts. Default 20."
        ),
    )
    parser_qualitysets.add_argument(
        "-f",
        "--qformat",
        dest="qformat",
        default="sanger",
        help=(
            'Quality format, options are "sanger" or "illumina". '
            "This will add an offset of 33 or 64 respectively to "
            'the qlimit value. Default "sanger".'
        ),
    )
    parser_qualitysets.add_argument(
        "-N",
        dest="n",
        type=int,
        default=20,
        help=(
            "Threshold to filter positions by the sum of read "
            "depth of the two samples. Default 20."
        ),
    )
    parser_genotype.add_argument(
        "--hom",
        dest="hom",
        type=float,
        default=0.9,
        help="Threshold to select homozygous positions. Default 0.9.",
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
    return parser.parse_args(argv)


def bam2seqz(subparsers, module_name, argv, log):
    log.log.info("Start %s" % module_name)
    args = bam2seqz_args(add_parser(subparsers, module_name), argv)
    if args.nproc > 1:
        if len(args.chr) < 2:
            msg = (
                "ERROR: The --chromosome parameter must include 2 or "
                "more chromosomes when used with the --parallel option"
            )
            log.log.error(msg)
            raise Exception(msg)
        else:
            partial_bam2seqz = partial(bam2seqz_main, args)
            pool = Pool(processes=args.nproc)
            pool.map_async(partial_bam2seqz, args.chr).get(9999999)
            pool.close()
            pool.join()

    else:
        bam2seqz_main(args, args.chr)


def bam2seqz_main(args, regions):
    if args.nproc > 1:
        if args.out == "-":
            raise Exception(
                (
                    "ERROR: The --output parameter must be different then "
                    "STDOUT when used with the --parallel option"
                )
            )
        else:
            prefix, extension = split_ext(args.out)
            file_output = "%s_%s%s" % (
                prefix,
                regions.replace(":", "_").replace("-", "_"),
                extension,
            )
            regions = [regions]
    else:
        file_output = args.out
    if args.pileup is False:
        program_bin = args.samtools
        bam = True
        if args.fasta is None:
            raise Exception(
                ("ERROR: The --fasta parameter is required " "when using BAM files")
            )
    else:
        program_bin = args.tabix
        bam = False

    qlimit = args.qlimit
    if args.qformat == "sanger":
        qlimit = qlimit + 33
    elif args.qformat == "illumina":
        qlimit = qlimit + 64
    else:
        raise Exception(
            ('Supported quality format are only "illumina" ' 'and "sanger"(default).')
        )

    bed = args.target_file
    if args.normal2:
        with open_bam_pileup(
            bam=bam,
            file=args.normal,
            fasta=args.fasta,
            program_bin=program_bin,
            regions=regions,
            bed=bed,
        ) as normal, open_bam_pileup(
            bam=bam,
            file=args.tumor,
            fasta=args.fasta,
            program_bin=program_bin,
            regions=regions,
            bed=bed,
        ) as tumor, open_bam_pileup(
            bam=bam,
            file=args.normal2,
            fasta=args.fasta,
            program_bin=program_bin,
            regions=regions,
            bed=bed,
        ) as alt_normal, xopen(
            args.gc, "rt"
        ) as gc_file, xopen(
            file_output, "wt", bgzip=True
        ) as out:
            mpup = zip_coordinates(split_coordinates(normal), split_coordinates(tumor))
            mpup_gc = zip_coordinates(mpup, Wiggle(gc_file))
            alt_mpup = zip_coordinates(mpup_gc, split_coordinates(alt_normal))
            write_seqz(
                input_zip=alt_mpup,
                output=out,
                depth_sum=args.n,
                qlimit=qlimit,
                hom_t=args.hom,
                het_t=args.het,
                het_f=args.het_f,
            )
    else:
        with open_bam_pileup(
            bam=bam,
            file=args.normal,
            fasta=args.fasta,
            program_bin=program_bin,
            regions=regions,
            bed=bed,
        ) as normal, open_bam_pileup(
            bam=bam,
            file=args.tumor,
            fasta=args.fasta,
            program_bin=program_bin,
            regions=regions,
            bed=bed
        ) as tumor, xopen(
            args.gc, "rt"
        ) as gc_file, xopen(
            file_output, "wt", bgzip=True
        ) as out:
            mpup = zip_coordinates(split_coordinates(normal), split_coordinates(tumor))
            mpup_gc = zip_coordinates(mpup, Wiggle(gc_file))
            write_seqz(
                input_zip=mpup_gc,
                output=out,
                depth_sum=args.n,
                qlimit=qlimit,
                hom_t=args.hom,
                het_t=args.het,
                het_f=args.het_f,
            )
    if file_output.endswith(".gz"):
        tabix_seqz(file_output, args.tabix)


def write_seqz(input_zip, output, depth_sum, qlimit, hom_t, het_t, het_f):
    """
    Wrap around the seqz output write to make the main function less verbose
    """
    header = seqz_header()
    output.write("\t".join(map(str, header)) + "\n")
    for line in input_zip:
        coordinate, data = line
        chromosome, position, null = coordinate
        seqz_line = do_seqz(
            data,
            depth_sum=depth_sum,
            qlimit=qlimit,
            hom_t=hom_t,
            het_t=het_t,
            het_f=het_f,
        )
        if seqz_line:
            output.write(
                "%s\t%i\t%s\n" % (chromosome, position, "\t".join(map(str, seqz_line)))
            )


def open_bam_pileup(bam, file, fasta, program_bin, regions, bed):
    """
    Wrap around a file object, switch the program used to open
    the file according to the bam flag and the regions options
    """
    if bam:
        return closing(
            bam_mpileup(
                bam=file, fasta=fasta, samtools_bin=program_bin, regions=regions, bed=bed
            )
        )
    else:
        if len(regions) == 0:
            return xopen(file, "rt")
        else:
            return closing(
                indexed_pileup(pileup=file, tabix_bin=program_bin, regions=regions, bed=bed)
            )
