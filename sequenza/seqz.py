from __future__ import division

try:
    import __pypy__
    from sequenza.pileup import acgt
except ImportError:
    try:
        from sequenza.c_pileup import acgt
    except ImportError:
        from sequenza.pileup import acgt


def seqz_header():
    header = [
        "chromosome",
        "position",
        "base.ref",
        "depth.normal",
        "depth.tumor",
        "depth.ratio",
        "Af",
        "Bf",
        "zygosity.normal",
        "GC.percent",
        "good.reads",
        "AB.normal",
        "AB.tumor",
        "tumor.strand",
    ]
    return header


def unpack_data(data):
    """
    Unpack normal, tumor and gc info from the
    specific touple structure and remove redundant information
    """
    if len(data) == 4:
        normal_line, tumor_line, gc, alt_line = data
        try:
            null, alt_depth, null = alt_line.split("\t", 2)
        except ValueError:
            alt_depth = 0
        alt_depth = int(alt_depth)

    else:
        alt_depth = None
        normal_line, tumor_line, gc = data

    ref, p1_str = normal_line.split("\t", 1)
    null, p2_str = tumor_line.split("\t", 1)

    p1_list = p1_str.split("\t")
    p2_list = p2_str.split("\t")
    ref = str(ref).upper()
    if alt_depth:
        p1_list[0] = alt_depth
    else:
        p1_list[0] = int(p1_list[0])
    p2_list[0] = int(p2_list[0])
    return (ref, p1_list, p2_list, alt_depth, gc)


def acgt_genotype(acgt_dict, freq_list, strand_list, hom_t, het_t, het_f, bases_list):
    """
    Return the alleles in the genotype
    """
    alleles = list()
    AB_max = sorted(freq_list, reverse=True)[:2]
    ref_idx = freq_list.index(AB_max[0])
    alleles.append(ref_idx)
    if AB_max[0] < hom_t and AB_max[1] >= het_t:
        alt_idx = [
            i
            for i in range(len(freq_list))
            if freq_list[i] == AB_max[1] and i != alleles[0]
        ][0]
        alt = bases_list[alt_idx]
        alt_fw_freq = strand_list[alt_idx] / acgt_dict[alt]
        if het_f < alt_fw_freq < (1 - het_f):
            alleles.append(alt_idx)
    return alleles


def acgt2seqz(
    normal,
    tumor,
    hom_t,
    het_t,
    het_f,
    normal_depth,
    tumor_depth,
    ref,
    gc,
    bases_list,
    het_only=False,
):
    normal_strand = normal["Z"]
    normal.pop("Z")
    sum_normal = sum(normal.values())
    if sum_normal > 0:
        tumor_strand = tumor["Z"]
        tumor.pop("Z")
        sum_tumor = sum(tumor.values())
        if sum_tumor > 0:
            normal_freq = [normal[x] / sum_normal for x in bases_list]
            tumor_freq = [tumor[x] / sum_tumor for x in bases_list]
            # Genotype normal
            alleles = acgt_genotype(
                normal, normal_freq, normal_strand, hom_t, het_t, het_f, bases_list
            )
            sum_tumor = sum(tumor.values())
            tumor_freq = [tumor[x] / sum_tumor for x in bases_list]
            if len(alleles) == 1 and not het_only:
                return parse_homoz(
                    tumor,
                    tumor_freq,
                    sum_tumor,
                    tumor_strand,
                    normal,
                    normal_depth,
                    tumor_depth,
                    alleles,
                    ref,
                    gc,
                    bases_list,
                )
            elif len(alleles) == 2:
                # Sort the genotype in the tumor
                if tumor_freq[alleles[0]] < tumor_freq[alleles[1]]:
                    alleles = list(reversed(alleles))
                return parse_heteroz(
                    tumor,
                    tumor_freq,
                    sum_tumor,
                    tumor_strand,
                    normal,
                    normal_depth,
                    tumor_depth,
                    alleles,
                    ref,
                    gc,
                    bases_list,
                )


def do_seqz(
    data, depth_sum=20, qlimit=53, hom_t=0.85, het_t=0.35, het_f=-0.1, het_only=False
):
    info = unpack_data(data)
    bases_list = ["A", "C", "G", "T"]
    if info[1][0] + info[2][0] >= depth_sum and len(info[1]) == len(info[2]):
        normal = acgt(info[1][1], info[1][2], info[1][0], info[0], qlimit)
        tumor = acgt(info[2][1], info[2][2], info[2][0], info[0], qlimit)
        return acgt2seqz(
            normal,
            tumor,
            hom_t,
            het_t,
            het_f,
            info[1][0],
            info[2][0],
            info[0],
            info[4],
            bases_list,
            het_only,
        )


def parse_homoz(
    tumor,
    tumor_freq,
    sum_tumor,
    tumor_strand,
    normal,
    normal_depth,
    tumor_depth,
    alleles,
    ref,
    gc,
    bases_list,
):
    i = alleles[0]
    no_zero_idx = [
        x
        for x in range(len(tumor_freq))
        if tumor_freq[x] > 0 and normal[bases_list[x]] == 0
    ]
    no_zero_bases = [
        "%s%s" % (bases_list[x], round(tumor_freq[x], 3)) for x in no_zero_idx if x != i
    ]
    if len(no_zero_bases) == 0:
        no_zero_bases = "."
        strands_bases = "0"
    else:
        no_zero_bases = ":".join(map(str, no_zero_bases))
        strands_bases = [
            "%s%s"
            % (bases_list[ll], round(tumor_strand[ll] / tumor[bases_list[ll]], 3))
            for ll in no_zero_idx
            if ll != i
        ]
        strands_bases = ":".join(map(str, strands_bases))
    homoz_tumor = tumor_freq[i]
    return [
        ref,
        normal_depth,
        tumor_depth,
        round(tumor_depth / normal_depth, 3),
        round(homoz_tumor, 3),
        0,
        "hom",
        gc,
        sum_tumor,
        bases_list[i],
        no_zero_bases,
        strands_bases,
    ]


def parse_heteroz(
    tumor,
    tumor_freq,
    sum_tumor,
    tumor_strand,
    normal,
    normal_depth,
    tumor_depth,
    alleles,
    ref,
    gc,
    bases_list,
):
    genotype = "".join([bases_list[idx] for idx in alleles])
    tumor_freqs = sorted([tumor_freq[i] for i in alleles], reverse=True)
    return [
        ref,
        normal_depth,
        tumor_depth,
        "%.3f" % (int(tumor_depth) / int(normal_depth)),
        round(tumor_freqs[0], 3),
        round(tumor_freqs[1], 3),
        "het",
        gc,
        sum_tumor,
        genotype,
        ".",
        "0",
    ]


def binned_seqz(seqz_file, window):
    def reset_window():
        return [None, 0, 0, 0, 0, 0, 0]

    def format_seqz_window(window):
        seqz_template = "%s\t%i\tN\t%i\t%i\t%s\t1.0\t0\thom\t%s\t%s\tN\t.\t0\n"
        return seqz_template % (
            window[0],
            window[1],
            int(round(window[2] / window[6], 0)),
            int(round(window[3] / window[6], 0)),
            round(window[4] / window[6], 3),
            int(round(window[5] / window[6], 0)),
            window[6],
        )

    def replace_gc(lines, gc):
        if lines == "":
            return lines
        else:
            lines = lines.rstrip().split("\n")
            new_gc = list()
            for line in lines:
                line = line.split("\t")
                line[9] = gc
                new_gc.append("\t".join(map(str, line)))
            return "%s\n" % "\n".join(new_gc)

    window_i = reset_window()
    buffer_line = ""
    for line in seqz_file:
        if line.startswith("chromosome"):
            yield line
        else:
            (
                chromosome,
                position,
                ref,
                depth_n,
                depth_t,
                depth_r,
                Af,
                Bf,
                genotype,
                gc,
                good_reads,
                alleles_n,
                alleles_t,
                strand_t,
            ) = line.split("\t")
            position = int(position)
            if genotype == "het" or alleles_t != ".":
                buffer_line += line
            if window_i[0] is None:
                window_i[0] = chromosome
                window_i[1] = int(position)
            if window_i[0] == chromosome:
                if position - window_i[1] <= window:
                    window_i[2] += int(depth_n)
                    window_i[3] += int(depth_t)
                    window_i[4] += float(depth_r)
                    window_i[5] += float(gc)
                    window_i[6] += 1
                else:
                    yield format_seqz_window(window_i) + replace_gc(
                        buffer_line, int(round(window_i[5] / window_i[6], 0))
                    )
                    buffer_line = ""
                    window_i = [
                        chromosome,
                        position,
                        int(depth_n),
                        int(depth_t),
                        float(depth_r),
                        float(gc),
                        1,
                    ]
            else:
                yield format_seqz_window(window_i) + replace_gc(
                    buffer_line, int(round(window_i[5] / window_i[6], 0))
                )
                window_i = [
                    chromosome,
                    position,
                    int(depth_n),
                    int(depth_t),
                    float(depth_r),
                    float(gc),
                    1,
                ]
                buffer_line = ""
    yield format_seqz_window(window_i) + replace_gc(
        buffer_line, int(round(window_i[5] / window_i[6], 0))
    )
