from sequenza.misc import xopen, countN
from sequenza.fasta import Fasta
from sequenza.wig import Wiggle


def add_parser(subparsers, module_name):
    return subparsers.add_parser(
        module_name,
        add_help=False,
        help=(
            "Given a fasta file and a window size it computes the GC "
            "percentage across the sequences, and returns a file in "
            "the UCSC wiggle format."
        ),
    )


def gc_wiggle_args(parser, argv):
    parser.add_argument(
        "-f",
        "--fasta",
        dest="fasta",
        required=True,
        help=("the fasta file. It can be a file " 'name or "-" to use STDIN'),
    )
    parser.add_argument(
        "-o", dest="out", type=str, default="-", help='Output file "-" for STDOUT'
    )
    parser.add_argument(
        "-w",
        dest="window",
        type=int,
        default=50,
        help=("The window size to calculate " "the GC-content percentage"),
    )
    return parser.parse_args(argv)


def gc_wiggle(subparsers, module_name, argv, log):
    log.log.info("Start %s" % module_name)
    args = gc_wiggle_args(add_parser(subparsers, module_name), argv)
    with xopen(args.fasta, "rt") as fa_file, xopen(args.out, "wt") as wg:
        stream_fasta = Fasta(fa_file, args.window)
        wiggle = Wiggle(wg)
        for seq in stream_fasta:
            if seq:
                nucleotides = seq[3]
                Ns = countN(nucleotides, "N")
                if Ns < 50:
                    if (seq[2] - seq[1] + 1) == args.window:
                        gc = int(
                            round(countN(nucleotides, "G"), 0)
                            + round(countN(nucleotides, "C"), 0)
                        )
                        wiggle.write((seq[0], seq[1], (seq[2] + 1), gc))
