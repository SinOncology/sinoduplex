
#include <htslib/sam.h>
#include "BamAlignment.h"


#include <sstream>



/// Initialize cigar_index_na, it is -1
const int BamAlignment::cigar_index_na = -1;

BamAlignment::BamAlignment()
  : isDuplicatedBamRecord(false)
{
  record = NULL;
  myAlignmentLength = 0;
}

/**
 * @brief constructor with  pointer to a bam record
 */

BamAlignment::BamAlignment(bam1_t *b)
  : isDuplicatedBamRecord(false)
{
  record = b;
  /// parse the cigar string and
  /// set the alignment length that this read spans in the reference
  parseCigarBinary();
}

/**
 * @brief constructor with a bam record and an indicator whether to duplicate the alignment or not
 */

BamAlignment::BamAlignment(bam1_t *b, bool duplicateBamRecord)
{
  if (duplicateBamRecord)
  {
    /// duplicate the bam record
    record = bam_dup1(b);
    isDuplicatedBamRecord = true;
  }
  else
  {
    /// just pointer to the bam record
    record = b;
    isDuplicatedBamRecord = false;
  }
  /// parse the cigar string and
  /// set the alignment length that this read spans in the reference
  parseCigarBinary();
  
}

/**
 * @brief De-constructor, clear the containers
 */
BamAlignment::~BamAlignment()
{
  
  if (isDuplicatedBamRecord)
  {
    bam_destroy1(record);
    record = NULL;
  }
  else
  {
    record = NULL;
  }
}

/**
 * @brief initiate from a bam1_t record
 * @param b
 */

void BamAlignment::readBamRecord(bam1_t *b)
{
  reset();
  record = b;
  isDuplicatedBamRecord = false;
  parseCigarBinary();
  
  
}

/**
 * @brief duplicate from a bam1_t record
 * @param b
 */
void BamAlignment::duplicateBamRecord(bam1_t *b)
{
  reset();
  /// duplicate the bam record
  record = bam_dup1(b);
  isDuplicatedBamRecord = true;
  parseCigarBinary();
}

void BamAlignment::reset()
{
  record = NULL;
  myAlignmentLength = 0;
  isDuplicatedBamRecord = false;
}

/**
 * @brief Return the query index associated with the specified reference offset
 * when the query starts at the specified reference position based on this cigar.
 * @param refOffset
 * @return the query index
 */
int BamAlignment::getQueryIndex(int32_t refOffset)
{
  return myCigarRoller.getQueryIndex(refOffset);
}

/**
 * @brief Return the query index associated with the specified reference position
 * when the query starts at the specified reference position based on this cigar.
 * @param refPosition The reference position
 * @param queryStartPos The start position of a read
 * @return the query index
 */
int BamAlignment::getQueryIndex(int32_t refPosition, int32_t queryStartPos)
{
  return myCigarRoller.getQueryIndex(refPosition, queryStartPos);
}

/**
 * @brief Return the reference position associated with the specified query index
 * when the query starts at the specified reference position based on this cigar.
 * @param queryIndex The query index in an alignment
 * @param queryStartPos The start position of an alignment
 * @return the reference position
 */
int32_t BamAlignment::getRefPosition(int32_t queryIndex, int32_t queryStartPos)
{
  return myCigarRoller.getRefPosition(queryIndex, queryStartPos);
}

/**
 * @brief Return the reference offset associated with the specified query index
 * when the query starts at the specified reference position based on this cigar.
 * @param queryIndex The query index in an alignment
 * @return the reference offset
 */
int32_t BamAlignment::getRefOffset(int32_t queryIndex)
{
  return myCigarRoller.getRefOffset(queryIndex);
}

/**
 * @brief Check the given reference position is  the head of this alignment
 * @param refPosition
 * @return true if this reference position is the head of this alignment, false if not
 */
bool BamAlignment::isHead(const int32_t refPosition)
{
  /// the start position of alignment record
  if (refPosition == record->core.pos)
  {
    return true;
  }
  else
  {
    return false;
  }
}

/**
 * @brief Check the given reference position is the tail of this alignment
 * @param refPosition
 * @return true if the position is the tail of this alignment, false if not
 */
bool BamAlignment::isTail(const int32_t refPosition)
{
  /// get the reference position (0-based) of the last base in this alignment
  /// it is  the start position of alignment + the expected reference bases spanned by this alignment -1
  int32_t lastBaseRefPosition = record->core.pos + myAlignmentLength - 1;
  if (refPosition == lastBaseRefPosition)
  {
    return true;
  }
  else
  {
    return false;
  }
  
}


/**
 * @brief Get 0-based inclusive leftmost position
 * @return 0-based inclusive leftmost position (start position)
 */
int32_t BamAlignment::get0BasedAlignmentStart()
{
  return record->core.pos;
}


/**
 * @brief Get 0-based inclusive rightmost position
 * @return 0-based inclusive rightmost position
 */
int32_t BamAlignment::get0BasedAlignmentEnd()
{
  /// if the alignment length is 0, return the start position
  if (myAlignmentLength == 0)
  {
    return record->core.pos;
  }
  return (record->core.pos + myAlignmentLength - 1);
  
}

/**
 *
 * @return The read name
 */
std::string BamAlignment::getReadName()
{
  return (std::string) bam_get_qname(record);
}

/**
 * @brief Get alignment length that a read spans
 * @return The alignment length of a read
 */
int BamAlignment::getAlignmentLength()
{
  return myAlignmentLength;
}

/**
 * @brief Get the read length
 * @return The read length
 */
int BamAlignment::getReadLength()
{
  return (int) record->core.l_qseq;
}

/**
 * @brief check whether there is a tag for barcode, if yes, parse it and store it
 *
 */
std::string BamAlignment::readBarcode()
{
  std::string barcode = "";
  uint8_t *bc = bam_aux_get(record, "BC"); //pointer to the barcode
  if (bc != NULL)
  {
    barcode = ((char *) bam_aux_get(record, "BC"));
    /// ignore "Z"
    barcode = barcode.substr(1);
  }
  return barcode;
}

std::string BamAlignment::originCigar()
{
  std::string originCigar = "";
  uint8_t *bc = bam_aux_get(record, "OC"); //pointer to the barcode
  if (bc != NULL)
  {
    originCigar = ((char *) bam_aux_get(record, "OC"));
    /// ignore "Z"
    originCigar = originCigar.substr(1);
  }
  return originCigar;
}

std::string BamAlignment::mateCigar()
{
  std::string mateCigar = "";
  uint8_t *bc = bam_aux_get(record, "MC"); //pointer to the barcode
  if (bc != NULL)
  {
    mateCigar = ((char *)bc);
    /// ignore "Z"
    mateCigar = mateCigar.substr(1);
  }
  return mateCigar;
}


/**
 * @brief get the cigar string
 * @return cigar string
 */
std::string BamAlignment::getCigarString()
{
  return myCigarRoller.getString();
}

/**
 * @brief get the expanded cigar string
 * @return expanded cigar string
 */
std::string BamAlignment::getExpandedCigarString()
{
  std::string ex_s;
  myCigarRoller.getExpandedString(ex_s);
  return ex_s;
}

/**
 * @brief Return the number of clips that are at the beginning of the cigar.
 * @return
 */
int BamAlignment::getNumBeginClips() const
{
  return myCigarRoller.getNumBeginClips();
}


/**
 * @brief Return the number of clips that are at the end of the cigar.
 * @return
 */
int BamAlignment::getNumEndClips() const
{
  return myCigarRoller.getNumEndClips();
}


char BamAlignment::getCigarCharOp(int32_t expandedCigarIndex)
{
  return myCigarRoller.getCigarCharOp(expandedCigarIndex);
}


char BamAlignment::getCigarCharOpFromQueryIndex(int32_t queryIndex)
{
  return myCigarRoller.getCigarCharOp(
    myCigarRoller.getExpandedCigarIndexFromQueryIndex(queryIndex));
}


char BamAlignment::getCigarCharOpFromRefOffset(int32_t refOffset)
{
  return myCigarRoller.getCigarCharOp(
    myCigarRoller.getExpandedCigarIndexFromRefOffset(refOffset));
}


char BamAlignment::getCigarCharOpFromRefPos(int32_t refPosition, int32_t queryStartPos)
{
  return myCigarRoller.getCigarCharOp(
    myCigarRoller.getExpandedCigarIndexFromRefPos(refPosition, queryStartPos));
}

/**
 * @brief Parse the cigar string and store the information in myCigarRoller
 */
void BamAlignment::parseCigarBinary()
{
  uint32_t *cigarPtr = bam_get_cigar(record);
  myCigarRoller.Set(cigarPtr, (uint16_t) record->core.n_cigar);
  //myCigarRoller.getCigarString(myCigar);
  myAlignmentLength = myCigarRoller.getExpectedReferenceBaseCount();
  myUnclippedStartOffset = myCigarRoller.getNumBeginClips();
  myUnclippedEndOffset = myCigarRoller.getNumEndClips();
  
}

/**
 * @brief Update the cigar array in bam record
 * @param cigar The new cigar array
 * @param n The number of cigar operations in the new cigar array
 */
void BamAlignment::updateCigarArray(uint32_t *cigar, int n)
{
  if (n != record->core.n_cigar)
  {
    int o = record->core.l_qname + record->core.n_cigar * 4;
    if (record->data_len + (n - record->core.n_cigar) * 4 > record->m_data)
    {
      record->m_data = record->data_len + (n - record->core.n_cigar) * 4;
      kroundup32(record->m_data);
      record->data = (uint8_t *) realloc(record->data, (size_t) record->m_data);
    }
    memmove(record->data + record->core.l_qname + n * 4, record->data + o, (size_t) (record->data_len - o));
    memcpy(record->data + record->core.l_qname, cigar, (size_t) (n * 4));
    record->data_len += (n - record->core.n_cigar) * 4;
    record->core.n_cigar = (uint32_t) n;
  }
  
  else
  {
    memcpy(record->data + record->core.l_qname, cigar, (size_t) (n * 4));
  }
  
}

/**
 * @brief write the bam record to a bam file
 * @param samFile
 */
void BamAlignment::writeToBam(samfile_t *samFile)
{
  samwrite(samFile, record);
}


/*
 * judge the molecular is covering the specified region
 */
bool
BamAlignment::pairIsCoveringRegion(std::string region_chr_name, const uint32_t region_start, const uint32_t region_end)
{
  std::string oc = "";
  int region_chr = -1;
  int chr = -1;
  int mate_chr = -1;
  if (region_chr_name == "chrX") region_chr = 23;
  else if (region_chr_name == "chrY") region_chr = 24;
  else
  {
    try
    {
      region_chr = std::stoi((region_chr_name.substr(region_chr_name.find("chr") + 3)).c_str()) - 1;
    }
    catch (std::exception e)
    {
      std::cerr << "no valid chrom numbre " << region_chr_name << std::endl;
      return -1;
    }
  }
  oc = originCigar();
  
  if (oc != "")
  {
    myCigarRoller.clear();
    myCigarRoller.Set(oc.c_str());
    myAlignmentLength = myCigarRoller.getExpectedReferenceBaseCount();
    myUnclippedStartOffset = myCigarRoller.getNumBeginClips();
    myUnclippedEndOffset = myCigarRoller.getNumEndClips();
  }
  chr = getReferenceID();
  mate_chr = getMateReferenceID();
  
  std::int32_t start = get0BasedAlignmentStart();
  std::int32_t end = get0BasedAlignmentEnd();
  std::int32_t mate_start = record->core.mpos;
  CigarRoller mate_cigar = CigarRoller(mateCigar().c_str());
  std::int32_t mate_end = mate_start + mate_cigar.getExpectedReferenceBaseCount();
  std::cout << "start: " << start << std::endl;
  std::cout << "end: " << end << std::endl;
  std::cout << "chr: " << chr << std::endl;
  std::cout << "region start: " << region_start << std::endl;
  std::cout << "region end: " << region_end << std::endl;
  std::cout << "region chr: " << region_chr << std::endl;
  std::cout << "mate start: " << mate_start << std::endl;
  std::cout << "mate end: " << mate_end << std::endl;
  std::cout << "mate chr: " << mate_chr << std::endl;
  
  
  bool result = false;
//  std::int32_t p1_start = -1;
//  std::int32_t p1_end = -1;
//  std::int32_t p2_start = -1;
//  std::int32_t p2_end = -1;
  //take the different chrom into consideration
  if (region_chr != chr && region_chr != mate_chr)
  {
    return false;
  }
  if (chr != mate_chr)
  {
    if (chr == region_chr)
    {
      result = (std::min(start, end) <= region_start) && (std::max(start, end) >= region_end);
      return result;
    }
    if (mate_chr == region_chr)
    {
      result = (std::min(mate_start, mate_end) <= region_start) && (std::max(mate_start, mate_end) >= region_end);
      return result;
      
    }
  }
  if (chr == mate_chr && mate_chr == region_chr)
  {
    result = ((std::min(start, end) <= region_start) && (std::max(start, end) >= region_end)) ||
             ((std::min(mate_start, mate_end) <= region_start) && (std::max(mate_start, mate_end) >= region_end));
    return result;
  }
  return result;
}



bool BamAlignment::readIsCoveringRegion(std::string region_chr_name, uint32_t region_start, uint32_t region_end)
{
  int region_chr = -1;
  int chr = -1;
  
  if (region_chr_name == "chrX") region_chr = 23;
  else if (region_chr_name == "chrY") region_chr = 24;
  else
  {
    try
    {
      region_chr = std::stoi((region_chr_name.substr(region_chr_name.find("chr") + 3)).c_str()) - 1;
    }
    catch (std::exception e)
    {
      std::cerr << "no valid chrom numbre " << region_chr_name << std::endl;
      return -1;
    }
  }
  chr = getReferenceID();
  bool result = false;
  if (region_chr != chr)
  {
    return false;
  }
  std::int32_t start = get0BasedAlignmentStart();
  std::int32_t end = get0BasedAlignmentEnd();
  result = (std::min(start, end) <= region_start) && (std::max(start, end) >= region_end);
  return result;
}

/*
 * return the covering len in a specified ref region in fact(when it cames D minus 1, when it cames I or M plus 1 )
 */

void BamAlignment::readCoveredLenInRefRepeatRegion(std::string region_chr, uint32_t region_start, uint32_t region_end,
                                                   int &form, int &length)
{
  int support_length = 0;
  std::string region_str = get_sequence_nib(region_chr, region_start, region_end);
  std::pair<int, std::string> Repeat = strRepeatEasy(region_str);
  if (Repeat.first == 0 || Repeat.second == "")
  {
//    length_form_pair = make_pair(NOT_REPEAT,0);
    form = NOT_REPEAT;
    length = 0;
    return;
  }
  
  // when the seq of target region is composed of repeat unit, let us compute the len of repeat units it supports
  
  if (getReadName() != region_chr)
  {
    form = ISOLATE;
    length = 0;
    return;
  }
  std::string readSeq = getReadSeq();
  std::string read_loci_repeat_seq = "";
  string readExpandedCigar, lociExpandedCigar;
  int lociStartReadIndex, lociEndReadIndex, match_count, del_count, insert_count;
  bool covering = false, isolate = false, partial = false;
  
  //TODO take bwa soft clipping into consideration
  //TODO take pair overlap into consideration
  std::pair<int32_t, int32_t> readSpanRefPos = make_pair(get0BasedAlignmentStart(), get0BasedAlignmentEnd());
  covering = (region_start >= readSpanRefPos.first) && (region_end <= readSpanRefPos.second);
  isolate = (region_start > readSpanRefPos.second) || (region_end < readSpanRefPos.first);
  partial = !isolate;
  if (covering)
  {
    readExpandedCigar = myCigarRoller.getExpandedString();
    lociStartReadIndex = myCigarRoller.getExpandedCigarIndexFromRefPos(region_start, readSpanRefPos.first);
    lociEndReadIndex = myCigarRoller.getExpandedCigarIndexFromRefPos(region_end, readSpanRefPos.first);
    lociExpandedCigar = readExpandedCigar.substr(lociStartReadIndex, lociEndReadIndex - lociStartReadIndex + 1);
    for (int i = 0; i < 2; ++i)
    {
    
    }
    match_count = countCharInStr(lociExpandedCigar, 'M');
    del_count = countCharInStr(lociExpandedCigar, 'D');
    insert_count = countCharInStr(lociExpandedCigar, 'I');
    form = EXACT_LENGTH;
    length = support_length;
    return;
  }
  if (partial)
  {
    form = SCOPE_LENGTH;
    length = support_length;
    return;
  }
  
  form = FAIL_TO_COMPUTE;
  length = 0;
  return;
  
}

std::string BamAlignment::getReadSeq()
{
//  static const char COMPLEMENT_BASE_MAP[16] = {0x0, 0x8, 0x4, 0xf, 0x2,0xf,0xf,0xf,0x1,0xf,0xf,0xf,0xf,0xf,0xf,0xf};
  std::string seq = "";
  seq = get_seq_str_bam(record);
  return seq;
}

long BamAlignment::getAlignmentStart()
{
  return get0BasedAlignmentStart() + 1;
}

long BamAlignment::getAlignmentEnd()
{
  return get0BasedAlignmentEnd() + 1;
}

std::string BamAlignment::getChrName()
{
  return chromID2ChrName(getReferenceID());
}

std::string BamAlignment::getSubExpandedCigarFromRefPos(int32_t refStartPos, int32_t refEndPos)
{
  int startCigarIndex, endCigarIndex;
  string subExpandedCigar = "";
  string expandedCigar = getExpandedCigarString();
  startCigarIndex = myCigarRoller.getExpandedCigarIndexFromRefPos(refStartPos, getAlignmentStart());
  endCigarIndex = myCigarRoller.getExpandedCigarIndexFromRefPos(refEndPos, getAlignmentStart());
  //int insert_base_num = 0;
  string tmp_str = "";
  Cigar cr;
  int j;
  if (startCigarIndex == -1 && endCigarIndex == -1)
  {
    return subExpandedCigar;
  }
  if (startCigarIndex == -1 && endCigarIndex >= 0)
  {
    subExpandedCigar = expandedCigar.substr(0, endCigarIndex + 1);
    return subExpandedCigar;
  }
  if (startCigarIndex >= 0 && endCigarIndex == -1)
  {
    subExpandedCigar = expandedCigar.substr(startCigarIndex, expandedCigar.size() - startCigarIndex + 1);
    if (subExpandedCigar.size() > (refEndPos - refStartPos + 1))
    {
      j = 0;
      for (int i = 0; i < refEndPos - refStartPos + 1;)
      {
        tmp_str += subExpandedCigar[j];
        if (subExpandedCigar[j] != 'I') i++;
        j++;
      }
      j = 0;
      subExpandedCigar = tmp_str;
    }
    return subExpandedCigar;
  }
  if (startCigarIndex >= 0 && endCigarIndex >= 0)
  {
    subExpandedCigar = expandedCigar.substr(startCigarIndex, endCigarIndex - startCigarIndex + 1);
    return subExpandedCigar;
  }
  return subExpandedCigar;
}

long BamAlignment::getAlignmentStartContainSoft()
{
  return get0BasedUnclippedStart() + 1;
}

long BamAlignment::getAlignmentEndContainSoft()
{
  return get0BasedUnclippedEnd() + 1;
}

std::string BamAlignment::getSubExpandedCigarFromRefPosContainSoft(int32_t refStartPos, int32_t refEndPos)
{
  
  
  return myCigarRoller.getCigarExpandedStringFromRefPosContainSoft(refStartPos, refEndPos,
                                                                   getAlignmentStartContainSoft());
}

void
BamAlignment::getInsertInfoFromRefPosContainSoft(int32_t refStartPos,
                                                 int32_t refEndPos, insert_info &target_insert)
{
  setInsertInfo();
  
  for (int i = 0; i < myInsertInfo.insertQueryIndex.size(); ++i)
  {
//    if
  }
}

void BamAlignment::setInsertInfo()
{
  myInsertInfo.insertBase.clear();
  myInsertInfo.insertLeftRefPos.clear();
  myInsertInfo.insertRightRefPos.clear();
  myInsertInfo.insertQueryIndex.clear();
  vector<int> queryIndexVec;
  queryIndexVec.clear();
  //int leftRefOffset=0;
  //int rightOffset=0;
  int leftQueryIndex =0;
  int rightQueryIndex = 0;
  myCigarRoller.getQueryIndexVecContainSoft(queryIndexVec);
  for (int i = 1; i <queryIndexVec.size() ; ++i)
  {
    if (myCigarRoller.getRefOffset(i) == -1){
      myInsertInfo.insertBase.push_back(getReadSeq()[i]);
      myInsertInfo.insertQueryIndex.push_back(i);
      myInsertInfo.insertLeftRefPos.push_back(leftQueryIndex);
      myInsertInfo.insertRightRefPos.push_back(rightQueryIndex);
    }else{
//      leftRefOffset
    }
  }
}

float BamAlignment::getMappingQual()
{
  return (float) record->core.qual;
}

/**
 * @brief read an alignment from bam file
 * @param fp
 * @param b
 * @param alignment
 * @return
 */



std::string BamAlignment::saTag()
{
  std::string  sa_tag;
  uint8_t *bc = bam_aux_get(record, "SA"); //pointer to the barcode
  if (bc != NULL)
  {
    sa_tag = ((char *)bc);
    /// ignore "Z"
    sa_tag = sa_tag.substr(1);
  }
  return sa_tag;
}


int readBamAlignment(samfile_t *fp, bam1_t *b, BamAlignment *alignment)
{
  int status = samread(fp, b);
  alignment->readBamRecord(b);
  return status;
}

