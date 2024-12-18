#ifndef __parsers_h_
#define __parsers_h_

void pileup_parse(int nucleot_list[4], int strand_list[4], char* pileup,
    char* quality, int depth, char reference, int qlimit,
    int noend, int nostart);

#endif
