#include "util_bam.h"

using namespace std;


int chrom2number(std::string chr)
{
  int i;
  if (chr == "chrX")
    i = 23;
  else if (chr == "chrY")
    i = 24;
  else
  {
    try
    {
      i = std::stoi((chr.substr(chr.find("chr") + 3)).c_str());
    }
    catch (exception e)
    {
      //std::cerr << "no valid chrom numbre " << chr << std::endl;
      return -1;
    }
  }
  return (i - 1);
}


/**
 * @brief combine chromsome ID and position into a 64 bit int
 * @param chromID
 * @param position
 * @return
 */

uint64_t combineChromPos(int32_t chromID, int32_t position)
{
  return (((uint64_t) chromID << 32) | (position & 0xFFFFFFFF));
}


string number2chrom(int num)
{
  // 0 is chr1, 21 is chr22, and so on
  if (num == 22) return "chrX";
  if (num == 23) return "chrY";
  return "chr" + std::to_string(num + 1);
}

/**
  * @brief Compute a text-based sequence string from a BAM structure
  * @param b The BAM structure pointer
  *
  * @returns A string containing the sequence
  */
std::string get_seq_str_bam(bam1_t *b)
{
  bam1_core_t *c = &b->core;
  std::stringstream line;
  
  uint8_t *ss = bam1_seq(b); //pointer to the sequence
  //Get the sequence
  line.str("");
  line.clear();
  for (int j = 0; j < c->l_qseq; ++j)
  {
    line << bam_nt16_rev_table[bam1_seqi(ss, j)];
  }
  
  return line.str();
}



/**
  * @brief Compute a text-based CIGAR string from a BAM structure
  * @param b The BAM structure pointer
  *
  * @returns A CIGAR string in text form
  */
std::string get_cigar_str_bam(bam1_t *b)
{
  
  uint32_t *cigar = bam1_cigar(b);
  int length = b->core.n_cigar;
  std::string str_cigar;
  std::stringstream line;
  line.str("");
  line.clear();
  int i = 0;
  // Go through all elements of the CIGAR string
  if (length > 0)
  {
    for (i = 0; i < length; i++)
    {
      // Collect length of constant CIGAR letter
      int cigar_pos_tmp = cigar[i] >> BAM_CIGAR_SHIFT;
      line << cigar_pos_tmp;
      // Collect current cigar letter
      line << BAM_CIGAR_STR[cigar[i] & BAM_CIGAR_MASK];
    }
    str_cigar = line.str();
  }
  else
    str_cigar = "*";
  
  return str_cigar;
}

std::string get_sequence_nib(std::string chrom, int32_t start_1based, int32_t end_1based)
{
  std::string seq="";
  //nib access
  nib nibObj;
  std::string build="hg19";
  //open genome access (nib)
  std::string fn=(std::string) nib_folder + "/" + build + "_" + chrom + ".nib";
  int stat=nibObj.open(fn);
  char base;
  for (int32_t i = start_1based ; i<= end_1based ; i++)
  {
    stat = nibObj.getBase(&base, i-1);
    seq +=base;
  }
  nibObj.close();
  return seq;
}


// refID is 0-based in bam file
std::string chromID2ChrName(int refID)
{
  std::string chromName = "";
  std::ostringstream os;
  if (refID == 23)
    chromName = "chrY";
  else if (refID == 22)
    chromName = "chrX";
  else if (refID >= 0 && refID < 22)
  {
    os << "chr" << refID + 1;
    chromName = os.str();
  }
  return (chromName);
}

std::string get_barcode(bam1_t *b)
{
  uint8_t *bc = bam_aux_get(b, "BC"); //pointer to the barcode
  /// Get barcodes of an alignment
  /// Complete line should be of the form: ZCGAGTAGTTGTATCCT (size 1 + 8 )
  std::string line((char *) bc);
  /// cut the first base Z
  return line.substr(1);
}


