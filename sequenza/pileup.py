def pileup_acgt(
    pileup, quality, depth, reference, qlimit=53, noend=False, nostart=False
):
    """
    Yet another version of the pileup parser. Used as a template
    for the C implementation, the old function still runs slightly
    faster, to my surprise...
    """
    characters = {
        "A": 0,
        "C": 1,
        "G": 2,
        "T": 3,
        "a": 4,
        "c": 5,
        "g": 6,
        "t": 7,
        "$": 8,
        "*": 9,
        "-": 10,
        "+": 11,
    }
    nucleot_list = [0, 0, 0, 0]
    strand_list = [0, 0, 0, 0]
    last_base = None
    index = len(pileup)
    n = 0
    q = 0
    while n < index:
        base = pileup[n]
        if base == "^":
            if not nostart:
                n += 2
                base = pileup[n]
            else:
                n += 3
                q += 1
                continue
        if base == ".":
            base = reference
        elif base == ",":
            base = reference.lower()
        base_index = characters[base]
        if base_index < 4:
            if ord(quality[q]) >= qlimit:
                last_base = base_index
                nucleot_list[base_index] += 1
                strand_list[base_index] += 1
            else:
                last_base = None
            n += 1
            q += 1
        elif base_index < 8:
            if ord(quality[q]) >= qlimit:
                last_base = base_index
                nucleot_list[base_index - 4] += 1
            else:
                last_base = None
            n += 1
            q += 1
        elif base_index == 8:
            if not noend:
                pass
            else:
                if last_base:
                    if last_base < 4:
                        nucleot_list[last_base] -= 1
                        strand_list[last_base] -= 1
                    else:
                        nucleot_list[last_base - 4] -= 1
            last_base = None
            n += 1
        elif base_index == 9:
            last_base = None
            n += 1
            q += 1
        elif base_index > 9:
            offset = n + 1
            while True:
                if pileup[offset].isdigit():
                    offset += 1
                else:
                    break
            step = int(pileup[n + 1 : offset])
            n = step + offset

    nucleot_dict = {
        "A": nucleot_list[0],
        "C": nucleot_list[1],
        "G": nucleot_list[2],
        "T": nucleot_list[3],
    }
    nucleot_dict["Z"] = [strand_list[0], strand_list[1], strand_list[2], strand_list[3]]
    return nucleot_dict


def acgt(pileup, quality, depth, reference, qlimit=53, noend=False, nostart=False):
    """
    Parse the mpileup format and return the occurrence of
    each nucleotides in the given positions.
    """
    nucleot_dict = {"A": 0, "C": 0, "G": 0, "T": 0}
    strand_dict = {"A": 0, "C": 0, "G": 0, "T": 0}
    n = 0
    block = {"seq": "", "length": 0}
    start = False
    del_ins = False
    l_del_ins = ""
    last_base = None
    ins_del_length = 0
    for base in pileup:
        if block["length"] == 0:
            if base == "$":
                if noend:
                    if last_base:
                        nucleot_dict[last_base.upper()] -= 1
                        if last_base.isupper():
                            strand_dict[last_base.upper()] -= 1
                    last_base = None
            elif base == "^":
                start = True
                block["length"] += 1
                block["seq"] = base
            elif base == "+" or base == "-":
                del_ins = True
                block["length"] += 1
                block["seq"] = base
            elif base == "." or base == ",":
                if ord(quality[n]) >= qlimit:
                    nucleot_dict[reference] += 1
                    if base == ".":
                        strand_dict[reference] += 1
                        last_base = reference
                    else:
                        last_base = reference.lower()
                else:
                    last_base = None
                n += 1
            elif base.upper() in nucleot_dict:
                if ord(quality[n]) >= qlimit:
                    nucleot_dict[base.upper()] += 1
                    if base.isupper():
                        strand_dict[base.upper()] += 1
                    last_base = base
                else:
                    last_base = None
                n += 1
            else:
                n += 1
        else:
            if start:
                block["length"] += 1
                block["seq"] += base
                if block["length"] == 3:
                    if not nostart:
                        if base == "." or base == ",":
                            if ord(quality[n]) >= qlimit:
                                nucleot_dict[reference] += 1
                                if base == ".":
                                    strand_dict[reference] += 1
                        elif base.upper() in nucleot_dict:
                            if ord(quality[n]) >= qlimit:
                                nucleot_dict[base.upper()] += 1
                                if base.isupper():
                                    strand_dict[base.upper()] += 1
                    block["length"] = 0
                    block["seq"] = ""
                    start = False
                    last_base = None
                    n += 1
            elif del_ins:
                if base.isdigit():
                    l_del_ins += base
                    block["seq"] += base
                    block["length"] += 1
                else:
                    ins_del_length = int(l_del_ins) + 1 + len(l_del_ins)
                    block["seq"] += base
                    block["length"] += 1
                    if block["length"] == ins_del_length:
                        block["length"] = 0
                        block["seq"] = ""
                        l_del_ins = ""
                        # ins_del = False
                        ins_del_length = 0

    nucleot_dict["Z"] = [
        strand_dict["A"],
        strand_dict["C"],
        strand_dict["G"],
        strand_dict["T"],
    ]
    return nucleot_dict
