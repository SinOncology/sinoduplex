#pragma once
#include "sam.h"
#include <string>
#include <vector>
#include <iostream>
#include <htslib/sam.h>
#include "CigarRoller.h"
#include "Cigar.h"

#include "util_bam.h"
#include "util_string.h"

#define NOT_REPEAT 0
#define ISOLATE 1
#define SCOPE_LENGTH 2
#define EXACT_LENGTH 3
#define FAIL_TO_COMPUTE 4

typedef struct
{
  vector<int32_t> insertLeftRefPos;
  vector<int32_t> insertRightRefPos;
  vector<int> insertQueryIndex;
  vector<char> insertBase;
} insert_info;

/// @brief An extend bam class to store more information about an alignment
class BamAlignment {
public:
  static const int cigar_index_na ;
  /// constructors
  BamAlignment();
  BamAlignment(bam1_t *b);
  BamAlignment(bam1_t *b,bool duplicateBamRecord);
  /// de-constructor
  ~BamAlignment ();
  /// initiate an alignment
  void readBamRecord(bam1_t *b);
  
  /// duplicate an alignment
  void duplicateBamRecord(bam1_t *b);
  
  void reset();
  bam1_t* getBamRecord(){return this->record; };
  
  /// Return the reference offset associated with the specified
  /// query index or INDEX_NA based on this cigar.
  int getRefOffset(int32_t queryIndex);
  
  /// Return the query index associated with the specified
  /// reference offset or INDEX_NA based on this cigar.
  int getQueryIndex(int32_t refOffset);
  
  /// Return the query index or INDEX_NA associated with the specified
  /// reference offset when the query starts at the specified reference
  /// position.
  int getQueryIndex(int32_t refPosition, int32_t queryStartPos);
  
  /// Return the reference position associated with the specified query index
  /// or INDEX_NA based on this cigar and the specified queryStartPos which
  /// is the leftmost mapping position of the first matching base in the
  /// query.
  int32_t getRefPosition(int32_t queryIndex, int32_t queryStartPos);
  
  /// check whether this alignment is head at the given reference position
  bool isHead(const int32_t refPosition);
  
  /// check whether this alignment is tail at the given reference position
  bool isTail(const int32_t refPosition);
  
  /// Returns the 0-based inclusive leftmost position of the alignment
  int32_t get0BasedAlignmentStart();
  
  /// Returns the 0-based inclusive rightmost position of the alignment
  int32_t get0BasedAlignmentEnd();
  
  /// get the read name
  std::string getReadName();
  
  /// get cigar string
  std::string getCigarString();
  
  /// get expanded cigar string
  std::string getExpandedCigarString();
  
  /// get the alignment length
  int getAlignmentLength() ;
  
  /// get the read length
  int getReadLength();
  
  /// get the barcode
  std::string readBarcode();
  std::string originCigar();
  std::string mateCigar();
  std::string saTag();
  
  /// Return the number of clips that are at the beginning of the cigar.
  int getNumBeginClips() const;
  
  /// Return the number of clips that are at the end of the cigar.
  int getNumEndClips() const;
  
  /// cigar operation related functions
  
  /// Return the character code of the cigar operator associated with the
  /// specified expanded CIGAR index.  '?' is returned for an out of range
  /// index.
  char getCigarCharOp(int32_t expandedCigarIndex);
  
  /// Return the character code of the cigar operator associated with
  /// the specified queryIndex.  '?' is returned for an out of range index.
  char getCigarCharOpFromQueryIndex(int32_t queryIndex);
  
  /// Return the character code of the cigar operator associated with
  /// the specified reference offset.  '?' is returned for an out of range offset.
  char getCigarCharOpFromRefOffset(int32_t refOffset);
  
  /// Return the character code of the cigar operator associated with
  /// the specified reference position.  '?' is returned for an out of
  /// range reference position.
  char getCigarCharOpFromRefPos(int32_t refPosition, int32_t queryStartPos);
  
  /// Parse the cigar string to calculate the cigar length and alignment end
  void parseCigarBinary();

  int32_t  getFlag(){return record->core.flag;};
  int32_t  getReferenceID() {return record->core.tid;}
  int32_t  getMateReferenceID(){return record->core.mtid;};
  
  int32_t  get0BasedPosition(){return record->core.pos;};
  int32_t  get0BasedMatePosition(){return record->core.mpos;};

  int32_t get0BasedUnclippedStart(){return (record->core.pos - myUnclippedStartOffset) ;};
  int32_t  get0BasedUnclippedEnd()
  {
    // myUnclippedEndOffset will be set by get0BasedAlignmentEnd if the
    // cigar has not yet been parsed, so no need to check it here.
    return (get0BasedAlignmentEnd() + myUnclippedEndOffset);
  }
  /// update the cigar array
  void updateCigarArray(uint32_t *cigar,int n);
  
  /// write the bam record to a bam file
   void writeToBam(samfile_t *samFile);
  
  void setFlag(uint16_t flag) {record->core.flag = flag;};
  
  //add by jin
  bool pairIsCoveringRegion(std::string region_chr, uint32_t region_start, uint32_t region_end);
  bool readIsCoveringRegion(std::string region_chr, uint32_t region_start, uint32_t region_end);
  void readCoveredLenInRefRepeatRegion(std::string region_chr, uint32_t region_start, uint32_t region_end,
                                       int &form, int &length);
  std::string getReadSeq();
  std::string getChrName();
  
  float getMappingQual();
  
  std::string getSubExpandedCigarFromRefPos(int32_t refStartPos, int32_t refEndPos);
  
  std::string getSubExpandedCigarFromRefPosContainSoft(int32_t refStartPos, int32_t refEndPos);
  long getAlignmentStart();
  
  long getAlignmentEnd();
  
  long getAlignmentStartContainSoft();
  
  long getAlignmentEndContainSoft();
  
  void getInsertInfoFromRefPosContainSoft(int32_t refStartPos,
                                          int32_t refEndPos, insert_info &target_insert);
  
  

protected:
  
  void setInsertInfo();
private:
  /// a pointer to the bam record containing information about an alignment in the bam file
  bam1_t *record;
  bool isDuplicatedBamRecord;
  /// The length of the alignment that the read spans in the reference genome.
  int myAlignmentLength;
  //std::string myCigar;
  CigarRoller myCigarRoller;
  // Unclipped alignment positions.
  int32_t myUnclippedStartOffset;
  int32_t myUnclippedEndOffset;
  insert_info myInsertInfo;
  
  
};


int readBamAlignment(samfile_t *fp,bam1_t *b,BamAlignment *alignment);



