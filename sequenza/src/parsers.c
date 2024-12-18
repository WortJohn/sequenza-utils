#include "parsers.h"

int get_index(char a)
{
    int i = 0;
    if (a == 'A') {
        i = 0;
    } else if (a == 'C') {
        i = 1;
    } else if (a == 'G') {
        i = 2;
    } else if (a == 'T') {
        i = 3;
    } else if (a == 'a') {
        i = 4;
    } else if (a == 'c') {
        i = 5;
    } else if (a == 'g') {
        i = 6;
    } else if (a == 't') {
        i = 7;
    } else if (a == '$') {
        i = 8;
    } else if (a == '*') {
        i = 9;
    } else if (a == 'N') {
        i = 10;
    } else if (a == 'n') {
        i = 11;
    } else if (a == '-') {
        i = 12;
    } else if (a == '+') {
        i = 13;
    }
    return i;
}

char get_lower(char nt) {
    if (nt == 'A') {
        nt = 'a';
    } else if (nt == 'C') {
        nt = 'c';
    } else if (nt == 'G') {
        nt = 'g';
    } else if (nt == 'T') {
        nt = 't';
    } else if (nt == 'N') {
        nt = 'n';
    }
   return nt;
}

int is_number (char a) {
    int i = 0;
    if (a == '0') {
        i = 1;
    } else if (a == '1') {
        i = 1;
    } else if (a == '2') {
        i = 1;
    } else if (a == '3') {
        i = 1;
    } else if (a == '4') {
        i = 1;
    } else if (a == '5') {
        i = 1;
    } else if (a == '6') {
        i = 1;
    } else if (a == '7') {
        i = 1;
    } else if (a == '8') {
        i = 1;
    } else if (a == '9') {
        i = 1;
    }
    return i;
}

void pileup_parse(int nucleot_list[4], int strand_list[4],
    char* pileup, char* quality, int depth, char reference,
    int qlimit, int noend, int nostart)
{
    char rev_reference;
    int n, q, last_base, offset, indel, step;
    n = 0;
    q = 0;
    last_base = -1;

   rev_reference = get_lower(reference);

    while(pileup[n] != '\0') {
        char base = pileup[n];
        if (base == '^') {
            if (nostart == 0) {
                n += 2;
                base = pileup[n];
            } else {
                n += 3;
                q += 1;
                continue;
            }
        }
        if (base == '.') {
            base = reference;
        } else if( base == ',') {
            base = rev_reference;
        }
        int index = get_index(base);
        int qual = quality[q];
        if (index < 4){
            if (qual >= qlimit) {
                last_base = index;
                nucleot_list[index] += 1;
                strand_list[index] += 1;
            } else {
                last_base = -1;
            }
            n += 1;
            q += 1;
        } else if (index < 8) {
            if (qual >= qlimit) {
                last_base = index;
                nucleot_list[index - 4] += 1;
            } else {
                last_base = -1;
            }
            n += 1;
            q += 1;
        } else if (index == 8) {
            if (noend == 1) {
                if (last_base > -1) {
                    if (last_base < 4) {
                        nucleot_list[last_base] -= 1;
                        strand_list[last_base]  -= 1;
                    } else {
                        nucleot_list[last_base - 4] -= 1;
                    }
                }
            }
            last_base = -1;
            n += 1;
        } else if (index >= 9 && index <= 11) {
            last_base = -1;
            n += 1;
            q += 1;
        } else if (index > 11) {
            offset = n + 1;
            indel = 0;
            for(;;) {
                if (is_number(pileup[offset])) {
                    if (indel == 0) {
                        step = pileup[offset] - '0';
                    } else {
                        step = (step * 10) + (int)(pileup[offset] - '0');
                    }
                    indel += 1;
                    offset += 1;
                } else {
                    break;
                }
            }
            n = step + offset;
            indel = 0;
        }
    }
}
