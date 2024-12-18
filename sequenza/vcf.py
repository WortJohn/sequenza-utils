from __future__ import division

import re
from sequenza.izip import zip_coordinates
from sequenza.seqz import acgt2seqz, seqz_header


def vcf_headline_content(line):
    """
    Try to get the string enclosed by "< ... >" in the VCF header
    """
    start = None
    end = None
    counter = 0
    for character in line:
        if start is None:
            if character == "<":
                start = counter + 1
        else:
            if character == ">":
                end = counter
        counter += 1
    if start and end:
        column = line[0 : start - 1].strip("##").strip("=")
        content = line[start:end]
        pattern = re.compile(r"""((?:[^,"']|"[^"]*"|'[^']*')+)""")
        content = pattern.split(content)[1::2]
        header_line = {}
        for item in content:
            key, value = item.split("=", 1)
            value = value.strip('"').split(":")
            if len(value) == 2:
                header_line[key] = [value[0], value[1].strip(" ").split("|")]
            else:
                header_line[key] = value[0]
        line_id = header_line.pop("ID", None)
        return (column, line_id, header_line)


def vcf_header(vcf_line, field="FORMAT", header=dict()):
    if vcf_line.startswith("##%s" % field):
        header_line = vcf_headline_content(vcf_line)
        if header_line:
            column = header_line[0]
            if column not in header.keys():
                header[column] = {}
                header[column][header_line[1]] = header_line[2]
            else:
                header[column][header_line[1]] = header_line[2]
    return header


def vfc2seqz(parsed_vcf, gc_wig, hom_t=0.85, het_t=0.35, het_f=-0.1, het_only=False):
    vcf_gc = zip_coordinates(parsed_vcf, gc_wig)
    header = seqz_header()
    yield header
    for line in vcf_gc:
        pos = line[0]
        vcf_line = line[1]
        gc = vcf_line[4]
        ref = vcf_line[0]
        alt = vcf_line[1].split(",")[0]
        normal = touple2acgt(ref, alt, vcf_line[2])
        tumor = touple2acgt(ref, alt, vcf_line[3])
        base_list = [ref, alt]
        if sum([normal[i] for i in base_list]) > 0 and vcf_line[2][0] > 0:
            seqz_line = acgt2seqz(
                normal,
                tumor,
                hom_t,
                het_t,
                het_f,
                vcf_line[2][0],
                vcf_line[3][0],
                ref,
                gc,
                [ref, alt],
                het_only,
            )
            if seqz_line:
                yield [pos[0], pos[1]] + seqz_line


def touple2acgt(ref, alt, item):
    item = list(item[1])
    if len(item) >= 4:
        a = item[0] + item[1]
        b = item[2] + item[3]
        try:
            z_a = item[0] / a
        except ZeroDivisionError:
            z_a = 0
        try:
            z_b = item[2] / b
        except ZeroDivisionError:
            z_b = 0
        Z = [z_a, z_b]
    else:
        a = item[0]
        b = item[1]
        Z = [0, 0]
    return {ref: a, alt: b, "Z": Z}


def vcf_parse(
    vcf_file,
    sample_order="n/t",
    field="FORMAT",
    depth=["DP", "DP"],
    alleles=["AD", "AD"],
    preset=None,
):
    """
    Parse the specified tags of a vcf file to retrieve total
    and per-allele depth information.
    """
    vcf_read = vcf_reader(vcf_file)
    header_list = dict()
    header = None
    tumor_col = None
    normal_col = None
    format_col = None
    for vcf_item in vcf_read:
        if "headerlist" in vcf_item:
            header_list = vcf_header(vcf_item["headerlist"], field, header_list)
        elif "header" in vcf_item:
            header = vcf_item["header"]
            last_col = len(header) - 1
            format_col = header.index(field)
            if sample_order == "n/t":
                tumor_col = last_col
                normal_col = tumor_col - 1
            elif sample_order == "t/n":
                normal_col = last_col
                tumor_col = normal_col - 1
        elif "line" in vcf_item:
            line_i = vcf_item["line"]
            chromosome = line_i[0]
            position = int(line_i[1])
            ref_alt = line_i[3:5]
            tumor_i = line_i[tumor_col]
            normal_i = line_i[normal_col]
            format_i = line_i[format_col]
            normal_i = split_format(
                format_i, normal_i, depth[0], alleles[0], ref_alt, ":", ",", preset
            )
            tumor_i = split_format(
                format_i, tumor_i, depth[1], alleles[1], ref_alt, ":", ",", preset
            )
            yield (
                (chromosome, position, position + 1),
                (line_i[3], line_i[4], normal_i, tumor_i),
            )


def split_format(
    format_str,
    format_content,
    depth,
    alleles,
    ref_alt,
    split_char1=":",
    split_char2=",",
    preset=None,
):
    format_str = format_str.split(split_char1)
    format_content = format_content.split(split_char1)
    if preset is None:
        allele_idx = format_str.index(alleles)
        alleles_val = list(map(int, format_content[allele_idx].split(split_char2)))
        try:
            depth_idx = format_str.index(depth)
            depth_val = int(format_content[depth_idx])
        except ValueError:
            depth_val = sum(alleles_val)
        return (depth_val, alleles_val)
    else:
        if preset == "caveman":
            return preset_caveman(format_str, format_content, ref_alt)
        elif preset == "dp4":
            return preset_dp4(format_str, format_content, ref_alt)
        elif preset == "strelka2_som":
            return preset_strelka2som(format_str, format_content, ref_alt)
        else:
            return split_format(
                format_str, format_content, depth, alleles, ref_alt, ":", ",", None
            )


def preset_caveman(format_str, format_content, ref_alt):
    idx_fa = format_str.index("FAZ")
    idx_fc = format_str.index("FCZ")
    idx_fg = format_str.index("FGZ")
    idx_ft = format_str.index("FTZ")
    idx_ra = format_str.index("RAZ")
    idx_rc = format_str.index("RCZ")
    idx_rg = format_str.index("RGZ")
    idx_rt = format_str.index("RTZ")
    fa, fc, fg, ft, ra, rc, rg, rt = [
        int(format_content[i])
        for i in [idx_fa, idx_fc, idx_fg, idx_ft, idx_ra, idx_rc, idx_rg, idx_rt]
    ]
    depth_val = sum([fa, fc, fg, ft, ra, rc, rg, rt])
    alleles_dict = {"A": [fa, ra], "C": [fc, rc], "G": [fg, rg], "T": [ft, rt]}
    alleles_val = alleles_dict[ref_alt[0]] + alleles_dict[ref_alt[1]]
    return (depth_val, alleles_val)

def preset_dp4(format_str, format_content, ref_alt):
    allele_idx = format_str.index("DP4")
    alleles_val = list(map(int, format_content[allele_idx].split(",")))
    depth_val = sum(alleles_val)
    return (depth_val, alleles_val)


def preset_strelka2som(format_str, format_content, ref_alt):
    tier_idx = 1
    idx_fa = format_str.index("AU")
    idx_fc = format_str.index("CU")
    idx_fg = format_str.index("GU")
    idx_ft = format_str.index("TU")
    fa, fc, fg, ft = [
        list(map(int, format_content[i].split(",")))
        for i in [idx_fa, idx_fc, idx_fg, idx_ft]
    ]
    depth_val = sum([fa[tier_idx], fc[tier_idx], fg[tier_idx], ft[tier_idx]])
    alleles_dict = {
        "A": fa[tier_idx],
        "C": fc[tier_idx],
        "G": fg[tier_idx],
        "T": ft[tier_idx],
    }
    alleles_val = [alleles_dict[ref_alt[0]], alleles_dict[ref_alt[1]]]
    return (depth_val, alleles_val)


def vcf_reader(vcf_file):
    for line in vcf_file:
        if line.startswith("##"):
            yield {"headerlist": line}
        elif line.startswith("#"):
            yield {"header": line.rstrip().split("\t")}
        else:
            yield {"line": line.rstrip().split("\t")}
