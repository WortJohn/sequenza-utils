from __future__ import division

from sequenza.izip import zip_coordinates
from sequenza.seqz import acgt2seqz, seqz_header


def copynumber4seqz(varscan_copynumber):
    header = next(varscan_copynumber).rstrip().split("\t")
    n_pos_col = header.index("num_positions")
    tumor_depth_col = header.index("tumor_depth")
    normal_depth_col = header.index("normal_depth")
    gc_col = header.index("gc_content")
    for line in varscan_copynumber:
        line = line.rstrip().split("\t")
        chromosome = line[0]
        start = int(line[1])
        end = int(line[2])
        tumor_depth = int(round(float(line[tumor_depth_col]), 0))
        normal_depth = int(round(float(line[normal_depth_col]), 0))
        good_reads = tumor_depth * int(line[n_pos_col])
        yield (
            (chromosome, start, end),
            (
                normal_depth,
                tumor_depth,
                round(tumor_depth / normal_depth, 3),
                good_reads,
                line[gc_col],
            ),
        )


def somatic2seqz(varscan_somatic, gc_wig, hom_t=0.9, het_t=0.25, het_f=-0.1):
    gc_varscan = zip_coordinates(varscan_somatic, gc_wig)
    header = seqz_header()
    yield header
    for line in gc_varscan:
        chromosome = line[0][0]
        position = line[0][1]
        varscan_data = line[1][0]
        gc = line[1][1]
        seqz_line = acgt2seqz(
            varscan_data["normal"],
            varscan_data["tumor"],
            hom_t,
            het_t,
            het_f,
            varscan_data["normal_depth"],
            varscan_data["tumor_depth"],
            varscan_data["ref"],
            gc,
            varscan_data["bases_list"],
        )
        yield [chromosome, position] + seqz_line


def parse_somatic(varscan_somatic):
    snp_header = next(varscan_somatic).rstrip().split("\t")

    ref_col = snp_header.index("ref")
    alt_col = snp_header.index("var")
    tumor_ref_plus_col = snp_header.index("tumor_reads1_plus")
    tumor_alt_plus_col = snp_header.index("tumor_reads2_plus")
    normal_ref_plus_col = snp_header.index("normal_reads1_plus")
    normal_alt_plus_col = snp_header.index("normal_reads2_plus")
    tumor_ref_minus_col = snp_header.index("tumor_reads1_minus")
    tumor_alt_minus_col = snp_header.index("tumor_reads2_minus")
    normal_ref_minus_col = snp_header.index("normal_reads1_minus")
    normal_alt_minus_col = snp_header.index("normal_reads2_minus")
    for line in varscan_somatic:
        line = line.rstrip().split("\t")
        chromosome = line[0]
        position = int(line[1])
        alt = line[alt_col]
        ref = line[ref_col]
        tumor_ref_plus = int(line[tumor_ref_plus_col])
        tumor_alt_plus = int(line[tumor_alt_plus_col])
        normal_ref_plus = int(line[normal_ref_plus_col])
        normal_alt_plus = int(line[normal_alt_plus_col])
        tumor_ref_minus = int(line[tumor_ref_minus_col])
        tumor_alt_minus = int(line[tumor_alt_minus_col])
        normal_ref_minus = int(line[normal_ref_minus_col])
        normal_alt_minus = int(line[normal_alt_minus_col])
        tumor_alt = tumor_alt_plus + tumor_alt_minus
        tumor_ref = tumor_ref_plus + tumor_ref_minus
        normal_alt = normal_alt_plus + normal_alt_minus
        normal_ref = normal_ref_plus + normal_ref_minus
        bases_list = [ref, alt]
        tumor_dict = {
            ref: tumor_ref,
            alt: tumor_alt,
            "Z": [
                try_division(tumor_ref_plus, tumor_ref),
                try_division(tumor_alt_plus, tumor_alt),
            ],
        }
        normal_dict = {
            ref: normal_ref,
            alt: normal_alt,
            "Z": [
                try_division(normal_ref_plus, normal_ref),
                try_division(normal_alt_plus, normal_alt),
            ],
        }
        yield (
            (chromosome, position, position + 1),
            (
                {
                    "normal": normal_dict,
                    "tumor": tumor_dict,
                    "normal_depth": normal_alt + normal_ref,
                    "tumor_depth": tumor_alt + tumor_ref,
                    "ref": ref,
                    "bases_list": bases_list,
                },
            ),
        )


def try_division(a, b):
    try:
        return a / b
    except ZeroDivisionError:
        return 0
