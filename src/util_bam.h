#pragma once

#include <string>
#include <sstream>
#include <algorithm>
#include <math.h>

#include <sam.h>
#include <htslib/kstring.h>
#include <htslib/hts.h>

#include "installdir.h"
#include "nibtools.h"
#include "util.h"


using namespace std;

int chrom2number(std::string chr);
uint64_t combineChromPos(int32_t chromID, int32_t position);
string number2chrom(int num);

std::string get_seq_str_bam(bam1_t *b);
std::string get_cigar_str_bam(bam1_t *b);

std::string get_sequence_nib(std::string chrom, int32_t start_1based, int32_t end_1based);

std::string chromID2ChrName(int refID);
std::string get_barcode(bam1_t *b);





