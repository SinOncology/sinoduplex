#pragma once

#include "SinoExecutable.h"
#include <vector>
#include <set>
#include <map>
#include <htslib/sam.h>
#include "htslib/kstring.h"

#include "BamAlignment.h"
#include "SamFlag.h"
#include "Logger.h"
#include "util_bam.h"
#include "util_string.h"
#include <math.h>
#ifdef __APPLE__
#include <memory>
#else
#include <tr1/memory>
#endif

class Duplex_Consensus : public SinoExecutable
{
public:
  static void printDuplex_ConsensusDescription(std::ostream &os);
  
  void printDescription(std::ostream &os);
  
  void printUsage(std::ostream &os);
  
  int execute(int argc, char **argv);
  
  virtual const char *getProgramName()
  {
    return ("Duplex_Consensus");
  }
  
  Duplex_Consensus() :
    myOneChrom(false), considerOne(false),
    lastCoordinate(-1), lastReference(-1)
  {
  }
  
  ~Duplex_Consensus();

private:
  
  struct PairInfo
  {
    int bc1;
    int bc2;
    int m1_size;
    int m2_size;
    uint32_t m1_flag;     ///< flag for the first mate
    uint32_t m2_flag;     ///< flag for the second mate
    string m_qname; ///< read name
    string m1_seq; ///< m1 sequence string
    string m2_seq; ///< m2 sequence string
    string m1_qual; ///< m1 quality string
    string m2_qual; ///< m2 quality string
    
    PairInfo()
    {
      bc1 = 0;
      bc2 = 0;
      m1_size = 0;
      m2_size = 0;
      m_qname = ""; ///< query name
      m1_flag = 0;     ///< flag for the first mate
      m2_flag = 0;     ///< flag for the second mate
      m1_seq = ""; ///< m1 sequence string
      m2_seq = ""; ///< m2 sequence string
      m1_qual = ""; ///< m1 quality string
      m2_qual = ""; ///< m2 quality string
    }
    
    ~PairInfo()
    {
    
    }
  };
  
  
  struct DupKey
  {
    int32_t reference;
    int32_t coordinate;
    bool orientation;
    
    DupKey()
      : reference(0), coordinate(0), orientation(false)
    {
    }
    
    DupKey(const DupKey &key)
      : reference(key.reference), coordinate(key.coordinate),
        orientation(key.orientation)
    {
    }
    
    DupKey(int32_t ref, int32_t coord, bool orient, uint32_t lib)
      : reference(ref), coordinate(coord), orientation(orient)
    {
    }
    
    inline void initKey(BamAlignment &record)
    {
      reference = record.getReferenceID();
      coordinate = record.get0BasedUnclippedStart();
      orientation = SamFlag::isReverse(record.getFlag());
      if (orientation)
      {
        /// Reverse, so get the unclipped end.
        coordinate = record.get0BasedUnclippedEnd();
      }
    }
    
    inline void initKey_bam(bam1_t *b)
    {
      BamAlignment record(b);
      reference = b->core.tid;
      coordinate = b->core.pos;
      orientation = SamFlag::isReverse(b->core.flag);
      if (orientation)
      {
        // Reverse, so get the unclipped end.
        coordinate = record.get0BasedUnclippedEnd();
      }
    }
    
    
    inline void copy(const DupKey &key)
    {
      reference = key.reference;
      coordinate = key.coordinate;
      orientation = key.orientation;
      
    }
    
    inline void cleanupKey(int32_t referenceID, int32_t coord)
    {
      reference = referenceID;
      coordinate = coord - CLIP_OFFSET;
      orientation = false;
      
    }
    
    inline bool operator<(const DupKey &key) const
    {
      if (reference == key.reference)
      {
        // same reference, so check the coordinate.
        if (coordinate == key.coordinate)
        {
          // Same coordinate, so check orientation.
          return (orientation < key.orientation);
        }
        // Same Ref, & different coordinates, so just check that.
        return (coordinate < key.coordinate);
      }
      // Different references, so less than is determined by
      // checking the reference only.
      return (reference < key.reference);
    }
    
    inline DupKey &operator=(const DupKey &key)
    {
      reference = key.reference;
      coordinate = key.coordinate;
      orientation = key.orientation;
      return (*this);
    }
    
    ~DupKey()
    {
    
    }
  };
  
  struct BarcodeKey
  {
    int bc1;
    int bc2;
    bool m1_isR1; //here only use m1_isR1
    //bool m2_isR2;
    /// when mate1 is R1 and mate2 is R2, then m1_isR1=1 and m2_isR2 =1
    /// when mate1 is R2 and mate2 is R1, then m1_isR1=0 and m2_isR2=0
    //std::string pair_orientation;
    
    inline BarcodeKey &operator=(const BarcodeKey &key)
    {
      bc1 = key.bc1;
      bc2 = key.bc2;
      m1_isR1 = key.m1_isR1;
      return (*this);
    }
    
    inline bool operator<(const BarcodeKey &key) const
    {
      
      if (bc1 == key.bc1)
      {
        if (bc2 == key.bc2)
        {
          
          return m1_isR1 < key.m1_isR1;
        }
        return bc2 < key.bc2;
      }
      
      return (bc1 < key.bc1);
    }
    
    BarcodeKey(int a, int b, uint32_t flag)
    {
      bc1 = a;
      bc2 = b;
      if (flag & BAM_FREAD1)
        m1_isR1 = true;
      else
        m1_isR1 = false;
    }
    
  };
  
  
  // Each read is assigned a key based on its referenceID, coordinate, orientation, and libraryID
  // This structure stores the two keys in a paired end read.
  struct PairedKey
  {
    DupKey key1;
    DupKey key2;
    
    PairedKey(DupKey k1, DupKey k2)
    {
      if (k2 < k1)
      {
        key1 = k2;
        key2 = k1;
      }
      else
      {
        key1 = k1;
        key2 = k2;
      }
    }
  };
  
  
  // Paired key comparison operator used for sorting paired end reads.
  struct PairedKeyComparator
  {
    inline bool operator()(const PairedKey &lhs, const PairedKey &rhs) const
    {
      if (lhs.key2 < rhs.key2) return true;
      if (rhs.key2 < lhs.key2) return false;
      return lhs.key1 < rhs.key1;
    }
  };
  
  struct MateData
  {
    int n_cigar;
    int l_qseq;
    int BC_index;
    uint32_t flag;
    std::string cigar_str;
    std::string readName;
    DupKey key;
    uint8_t *seq;
    uint8_t *qual;
    
    
    MateData()
    {
      n_cigar = 0;
      seq = NULL;
      qual = NULL;
      
    }
    
    MateData(bam1_t *b, DupKey key1, int BC_index1, bool stored_readname)
    {
      key = key1;
      if (stored_readname)
        readName = (std::string) bam_get_qname(b);
      else
        readName = "";
      seq = (uint8_t *) malloc((b->core.l_qseq + 1) / 2);
      memcpy(seq, bam_get_seq(b), (b->core.l_qseq + 1) / 2);
      l_qseq = b->core.l_qseq;
      qual = (uint8_t *) malloc((b->core.l_qseq));
      memcpy(qual, bam_get_qual(b), (b->core.l_qseq));
      
      BC_index = BC_index1;
      flag = b->core.flag;
      n_cigar = b->core.n_cigar;
      cigar_str = get_cigar_str_bam(b);
    }
    
    ~MateData()
    {
      free(seq);
      free(qual);
      
    }
    
  };
  
  struct PairedData
  {
    shared_ptr<MateData> m1;
    shared_ptr<MateData> m2;
    
    PairedData(shared_ptr<MateData> p1, shared_ptr<MateData> p2)
    {
      m1 = p1;
      m2 = p2;
    }
    
    ~PairedData()
    {
    
    }
    
  };
  
  std::vector<std::string> barcodes_list_sino;
  std::vector<std::string> barcodes_list_idt;
  
  
  /// Map for storing reads until the mate is found.
  typedef std::multimap<uint64_t, shared_ptr<MateData> > MateMap;
  
  typedef std::pair<MateMap::iterator, bool> MateMapInsertReturn;
 // std::tr1::shared_ptr<MateMap> myMateMap(new(MateMap));
  
  MateMap myMateMap;
  //MateMap * myMateMap =new(MateMap);
  
  /// A multi-map from the key of a paired read to its read data (represent one family)
  typedef std::multimap<PairedKey, PairedData, PairedKeyComparator> PairedMap;
  typedef PairedMap::iterator PairedMapInsertReturn;
  typedef std::pair<PairedMap::iterator, PairedMap::iterator> FamilyIterator;
  
  PairedMap myPairedMap;
  
  typedef std::multimap<BarcodeKey, PairedData> BarcodePairedMap;
  
  std::map<int, int> family2size;
  std::map<int, int> DS_N_size;
  uint32_t reads_correct_barcodes = 0;
  int32_t reads_with_wrong_barcodes = 0;
  uint32_t read_pairs_sscs = 0;
  uint32_t dscs_reads = 0;
  uint32_t sscs_reads = 0;
  uint32_t reads_2_fastq = 0;
  uint32_t sscs_families = 0;
  uint32_t total_families = 0;
  uint32_t n_ignored_family = 0;
  //uint32_t n_ignored_family_cigar = 0;
  uint32_t n_pass_mate = 0;
  bool myOneChrom;
  bool considerOne;
  bool onlySingleStrand;
  int lastCoordinate;
  int lastReference;
  
  int myMinQual;
  
  float N_percent_cutoff = 0.5;
  uint32_t N_SS_family = 0;
  uint32_t N_DS_family = 0;
  FILE *R1_Fastq_Ptr;
  FILE *R2_Fastq_Ptr;
  
  
  static const int DEFAULT_MIN_QUAL;
  static const uint32_t CLIP_OFFSET;
  
  const int familysize_cutoff = 1;
  //const int cigar_count_cutoff = 1;
  const float abundant_nucleotide_cutoff = 0.75;
  
  static const int frscore_sscs = 30;
  static const int frscore_N = 10;
  static const int frscore_dscs = 60;;
  static const int frscore_N_dscs = 10;
  
  static const int phred_score = 33;
  const  int base5[5] = {1, 2, 4, 8, 15};
  
  string Nucl[5] = {"A", "C", "G", "T", "N"};
  
  // Once record is read, look back at previous reads and determine
  // if any no longer need to be kept for duplicate checking.
  // Call with NULL to cleanup all records.
  void cleanupPriorReads(bam1_t *b);
  
  // Determines if the current position has changed when we read record
  bool hasPositionChanged(bam1_t *b);
  
  // When a record is read, check if it is a duplicate or
  // store for future checking.
  void checkDuplex(bam1_t *b);
  
  bool isDuplex(PairInfo &p1, PairInfo &p2);
  
  
  
  
  bool form_ss_consensus(std::vector<uint8_t *> read_seq_m1,
                         std::vector<uint8_t *> read_seq_m2,
                         std::vector<uint8_t *> read_qual_m1,
                         std::vector<uint8_t *> read_qual_m2,
                         PairInfo &sscs, int read_length);
  
  void form_ss_consensus_phred(std::vector<uint8_t *> read_seq_m1,
                         std::vector<uint8_t *> read_seq_m2,
                         std::vector<uint8_t *> read_qual_m1,
                         std::vector<uint8_t *> read_qual_m2,
                         PairInfo &sscs, int read_length);
  
  void handleFamilyPair(PairedMap::iterator familyIterator);
  
  void cleanFamilyPair(PairedMap::iterator familyIterator);
  
  void fastq_out(FILE *f, PairInfo &read_pair, int mate, bool append12);
  
  bool generate_ds_consensus_seq(PairInfo &top, PairInfo &bottom, PairInfo &dscs);
  
  bool generate_ds_consensus_seq_phred(PairInfo &top, PairInfo &bottom, PairInfo &dscs);
  
  double_t get_error_pr( int givenBase,int base, int baseQ);
  
};

