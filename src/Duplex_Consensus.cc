#include <iostream>
#include <algorithm>
#include <unistd.h>
#include <getopt.h>
#include "Duplex_Consensus.h"
#include "installdir.h"



const int Duplex_Consensus::DEFAULT_MIN_QUAL = 15;
const uint32_t Duplex_Consensus::CLIP_OFFSET = 1000;

Duplex_Consensus::~Duplex_Consensus()
{
  /// clean up the maps.
  /// Free any paired records.
  myPairedMap.clear();
  /// Free any entries in the mate map.
  myMateMap.clear();
}

void Duplex_Consensus::printDuplex_ConsensusDescription(std::ostream &os)
{
  os << "duplex - make duplex consensus reads\n";
}


void Duplex_Consensus::printDescription(std::ostream &os)
{
  printDuplex_ConsensusDescription(os);
}


void Duplex_Consensus::printUsage(std::ostream &os)
{
  os << "Usage: ./sinotools duplex -i <InputBamFile> -o <OutputPrefix> [-c <chrom>][-s]"
     << "[-l <LogFile>] [-q <base quality>] [-n <N cutoff>] [-single] [-one]";
  
  os << std::endl << std::endl;
  os << "Required parameters :" << std::endl;
  os << "\t-i <InputBamFile> : Input BAM file name (must be sorted)" << std::endl;
  os << "\t-o <OutputPrefix> : Output file prefix" << std::endl;
  os << "\t-c <Chromosome>   : chromosome to read" << std::endl;
  os << "\t-single           : only do single strand consensus  " << std::endl;
  os << "\t-one              : also consider families of size 1" << std::endl;
  
}

int Duplex_Consensus::execute(int argc, char **argv)
{
  std::string inFile, outFile, logFile;
  bool verboseFlag = true;
  onlySingleStrand = false;
  myOneChrom = false;
  considerOne = false;
  myMinQual = DEFAULT_MIN_QUAL;
  int longindex;
  int opt;
  std::string chrom = "";
  
  // option definition
  static struct option longopts[] = {
    {"help", 0, 0, 'h'},
    {"help", 0, 0, '?'},
    {"i",    1, 0, 1},
    {"o",    1, 0, 2},
    {"n",    1, 0, 3},
    {"q",    1, 0, 4},
    {"c",    1, 0, 5},
    {"l",    1, 0, 6},
    {"single", 0, 0, 7},
    {"one",  0, 0, 8},
    {0,      0, 0, 0}
  };
  
  optind = 0;
  // parse command line arguments
  while ((opt = getopt_long_only(argc, argv, "h?", longopts, &longindex)) != -1)
  {
    switch (opt)
    {
    case 'h':
      printUsage(std::cerr);
      exit(1);
    case '?':
      printUsage(std::cerr);
      exit(1);
    case 1:
      inFile = (std::string) optarg;
      break;
    case 2:
      outFile = (std::string) optarg;
      break;
    case 3:
      N_percent_cutoff = atof(optarg);
      break;
    case 4:
      myMinQual = atoi(optarg);
      break;
    case 5:
      chrom = (std::string) optarg;
      myOneChrom = true;
      break;
    case 6:
      logFile = (std::string) optarg;
      break;
    case 7:
      onlySingleStrand = true;
      break;
    case 8:
      considerOne = true;
      break;
    default:
      std::cerr << "Error: cannot parse arguments.\n";
      printUsage(std::cerr);
      exit(1);
    }
  }
  if (inFile.size() == 0)
  {
    printUsage(std::cerr);
    std::cerr << "Specify an input bam file" << std::endl;
    return -1;
  }
  
  if (outFile.size() == 0)
  {
    printUsage(std::cerr);
    
    std::cerr << "Specify an output file" << std::endl;
    return -1;
  }
  
  if (logFile.size() == 0)
  {
    logFile = outFile + ".log";
  }
  
  dscs_reads = 0;
  Logger::gLogger = new Logger(logFile.c_str(), verboseFlag);
  R1_Fastq_Ptr = fopen((outFile + "_R1.fastq").c_str(), "w");
  R2_Fastq_Ptr = fopen((outFile + "_R2.fastq").c_str(), "w");
  // Array containing the list of barcodes
  read_barcode_list(barcodes_file_sino, barcodes_list_sino);
  read_barcode_list(barcodes_file_idt, barcodes_list_idt);
  
  samfile_t *fp = NULL;
  // read the bam file
  fp = samopen(inFile.c_str(), "rb", 0);
  if (fp == NULL)
  {
    std::cerr << "Error: can not open input bam-file.\n";
    std::cerr << inFile << endl;
    return (-1);
  }
  int chromID = chrom2number(chrom);
  
  bam_index_t *idx_bam = bam_index_load(inFile.c_str());
  
  if (idx_bam == NULL)
  {
    std::cerr << "Error: can not open bai of input bam-file.\n";
    std::cerr << inFile << endl;
    std::cerr << "Please index it first \n";
    return (-1);
  }
  /// chromosome start and end position
  int beg_t, end_t, tid_t;
  bam_parse_region(fp->header, chrom.c_str(), &tid_t, &beg_t, &end_t);
  //cout<<"chromID "<<tid_t<<"\t"<<beg_t<<"\t"<<end_t<<endl;
  /// Create bam iterator
  bam_iter_t iter_bam = bam_iter_query(idx_bam, tid_t, beg_t, end_t);
  cout<<"only do single strand consensus\t"<<onlySingleStrand<<endl;
  cout << "reading bam file  " << chrom << endl;
  cout << inFile << endl;
  
  lastReference = -1;
  lastCoordinate = -1;
  
  /// for keeping some basic statistics
  uint32_t total_reads = 0;
  uint32_t total_mapped_reads = 0;
  uint32_t valid_reads = 0;
  
  int32_t complete_overlap_mapped_reads = 0;
  int32_t numLog = 5000000;
  
  bam1_t *b = bam_init1();
  int32_t t_start, t_end, dt;
  t_start = time(NULL);
  
  uint32_t rflag_filter = BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FDUP;
  /// Now we start reading records
  //while (samread(fp, b) >= 0)
  while (bam_iter_read(fp->x.bam, iter_bam, b) >= 0)
  {
    /// here is only for chr0, to parse those read pairs mapped in different
    /// chromosomes
    if (myOneChrom && chromID == -1 && b->core.tid == b->core.mtid)
      continue;
    const bool m1_is_mapped = !(b->core.flag & BAM_FUNMAP);
    const bool m2_is_mapped = !(b->core.flag & BAM_FMUNMAP);
    
    if (!(b->core.flag & BAM_FSECONDARY))
    {
      total_reads++;
      if (m1_is_mapped)
        total_mapped_reads++;
    }
    if (total_reads % numLog == 0)
    {
      t_end = time(NULL);
      dt = t_end - t_start;
      cerr << "processed " << total_reads << " reads, ";
      cerr << "used time: " << dt / 60 << " min " << dt - 60 * (dt / 60) << " sec  \n";
    }
    /// ignore those unmapped reads and secondary alignments
    if (b->core.flag & rflag_filter)
      continue;
    /// ignore those mapped reads with mapping quality =0
    if (m1_is_mapped && b->core.qual == 0)
      continue;
    
    /// both reads mapped and  complete overlap
    if (m1_is_mapped && m2_is_mapped && b->core.pos == b->core.mpos &&
        b->core.tid == b->core.mtid)
    {
      complete_overlap_mapped_reads++;
    }
    
    /// those read pairs are both mapped and both with  mapq >0 go to sscs family
    if (m1_is_mapped && m2_is_mapped && b->core.qual > 0)
    {
      valid_reads++;
      
      /// if we have moved to a new position, look back
      /// at previous read pairs for duplex
      if (hasPositionChanged(b))
      {
        cleanupPriorReads(b);
      }
      /// put the record in the appropriate maps:
      /// mate1 in the myMateMap first
      /// paired reads go in myPairedMap when reach mate2
      checkDuplex(b);
    }
  }
  
  /// we're finished reading the bam file, so clean up the map
  cleanupPriorReads(NULL);
  cout << "Left mateMap\t" << myMateMap.size() << "\t" << myPairedMap.size() << endl;
  myMateMap.erase(myMateMap.begin(), myMateMap.end());
  
  //uint32_t total_final_reads = sscs_reads + dscs_reads;
  
  ofstream out;
  out.open((outFile + ".log").c_str());
  if (!out.is_open())
  {
    cerr << "Error: cannot output file" << outFile << ".log" << endl;
    exit(1);
  }
  
  out << "total reads\t" << total_reads << endl;
  out << "total mapped reads\t" << total_mapped_reads << "\tpercentage:\t" << total_mapped_reads * 100.0 / total_reads
      << endl;
  out << "complete overlap mapped reads\t" << complete_overlap_mapped_reads << "\tpercentage:\t"
      << complete_overlap_mapped_reads * 100.0 / total_reads << endl;
  out << "mapped reads go to barcode checking\t" << valid_reads << "\tpercentage:\t"
      << valid_reads * 100.0 / total_reads << endl;
  out << "pass mates\t" << n_pass_mate << "\tpercentage:\t"
      << n_pass_mate * 100.0 / total_reads << endl;
  
  out << "mapped reads with correct barcodes\t" << reads_correct_barcodes << "\tpercentage:\t"
      << reads_correct_barcodes * 100.0 / total_reads << endl;
  out << "mapped reads with wrong barcodes\t" << reads_with_wrong_barcodes << "\tpercentage:\t"
      << reads_with_wrong_barcodes * 100.0 / total_reads << endl;
  
  out << "mapped read pairs go to sscs family\t" << read_pairs_sscs * 2 << "\tpercentage:\t"
      << read_pairs_sscs * 200.0 / total_reads << endl;
  
  out << "total families\t" << total_families << endl;
  out << "sscs family ratio\t" << read_pairs_sscs * 1.0 / total_families << endl;
  
  out << "sscs families (size <= " << familysize_cutoff << ")\t" << n_ignored_family
      << " percentage:\t" << n_ignored_family * 100.0 / total_families << endl;
  out<<"these families will be considered\t"<<considerOne<<endl;
  out << "N_DS_family (>= " << N_percent_cutoff << ")\t" << N_DS_family << "\t"
      << " percentage:\t" << N_DS_family * 100.0 / (dscs_reads) << endl;;
  
  out << "sscs read pairs\t" << sscs_reads << "\tpercentage:\t" << sscs_reads * 100.0 / total_families << endl;
  out << "dscs read pairs\t" << dscs_reads << "\tpercentage:\t" << dscs_reads * 100.0 / total_families << endl;
  //out << "ratio of (sscs+dscs)/dscs\t" << total_families * 1.0 / dscs_reads << endl;
  //out << "percentage of (sscs*2)/total_sscs\t" << sscs_reads * 200.0 / total_families << endl;
  out << "percentage of (dscs*2)/total_sscs\t" << dscs_reads * 200.0 / total_families << endl;
  out << "total dscs+sscs read pairs left to final fastq files\t" << reads_2_fastq << "\t" << reads_2_fastq * 2
      << endl;
  out << "duplex consensus ratio\t" << read_pairs_sscs * 1.0 / reads_2_fastq << endl;
  
  
  out.close();
  
  cout << "total reads\t" << total_reads << endl;
  cout << "total mapped reads\t" << total_mapped_reads << "\tpercentage:\t" << total_mapped_reads * 100.0 / total_reads
       << endl;
  cout << "complete overlap mapped reads \t" << complete_overlap_mapped_reads << "\tpercentage:\t"
       << complete_overlap_mapped_reads * 100.0 / total_reads << endl;
  cout << "mapped reads go to barcode checking\t" << valid_reads << "\tpercentage:\t"
       << valid_reads * 100.0 / total_reads << endl;
  cout << "pass mates\t" << n_pass_mate << "\tpercentage:\t"
       << n_pass_mate * 100.0 / total_reads << endl;
  cout << "mapped reads with correct barcodes\t" << reads_correct_barcodes << "\tpercentage:\t"
       << reads_correct_barcodes * 100.0 / total_reads << endl;
  cout << "mapped reads with wrong barcodes\t" << reads_with_wrong_barcodes << "\tpercentage:\t"
       << reads_with_wrong_barcodes * 100.0 / total_reads << endl;
  
  
  cout << "mapped read pairs go to sscs family\t" << read_pairs_sscs * 2 << "\tpercentage:\t"
       << read_pairs_sscs * 200.0 / total_reads << endl;
  
  cout << "total families\t" << total_families << endl;
  cout << "sscs family ratio\t" << read_pairs_sscs * 1.0 / total_families << endl;
  
  cout << "sscs families (size <= " << familysize_cutoff << ")\t" << n_ignored_family
       << " percentage:\t" << n_ignored_family * 100.0 / total_families << endl;
  cout<<"these families will be considered\t"<<considerOne<<endl;
  
  cout << "N_DS_family (>= " << N_percent_cutoff << ")\t" << N_DS_family << "\t"
       << " percentage:\t" << N_DS_family * 100.0 / (dscs_reads) << endl;
  
  cout << "sscs read pairs\t" << sscs_reads << "\tpercentage:\t" << sscs_reads * 100.0 / total_families << endl;
  cout << "dscs read pairs\t" << dscs_reads << "\tpercentage:\t" << dscs_reads * 100.0 / total_families << endl;
 //cout << "ratio of (sscs+dscs)/dscs\t" << total_families * 1.0 / dscs_reads << endl;
  //cout << "percentage of (sscs*2)/total_sscs\t" << sscs_reads * 200.0 / total_families << endl;
  
  cout << "percentage of (dscs*2)/total_sscs\t" << dscs_reads * 200.0 / total_families << endl;
  cout << "total dscs+sscs read pairs left to final fastq files\t" << reads_2_fastq << "\t" << reads_2_fastq * 2
       << endl;
  cout << "duplex consensus ratio\t" << read_pairs_sscs * 1.0 / reads_2_fastq << endl;
  
  
  ofstream out_count;
  out_count.open((outFile + "_families_count.txt").c_str());
  if (!out_count.is_open())
  {
    cerr << "Error: cannot output file" << outFile << "_families_count.txt" << endl;
    exit(1);
  }
  for (map<int, int>::iterator it_s_ontarget = family2size.begin();
    it_s_ontarget != family2size.end(); it_s_ontarget++)
  {
    out_count << it_s_ontarget->first << "\t" << it_s_ontarget->second
              << "\t" << it_s_ontarget->first * it_s_ontarget->second << endl;
  }
  out_count.close();
  //
  out_count.open((outFile + "_N_count.txt").c_str());
  if (!out_count.is_open())
  {
    cerr << "Error: cannot output file" << outFile << "_N_count.txt" << endl;
    exit(1);
  }
  uint32_t total_N_reads = 0;
  for (map<int, int>::iterator it_s_ontarget = DS_N_size.begin();
    it_s_ontarget != DS_N_size.end(); it_s_ontarget++)
  {
    total_N_reads += it_s_ontarget->second;
  }
  
  for (map<int, int>::iterator it_s_ontarget = DS_N_size.begin();
    it_s_ontarget != DS_N_size.end(); it_s_ontarget++)
  {
    out_count << it_s_ontarget->first << "\t" << it_s_ontarget->second << "\t"
              << it_s_ontarget->second * 100.0 / total_N_reads << endl;
  }
  
  out_count.close();
  samclose(fp);
  bam_destroy1(b);
  fclose(R1_Fastq_Ptr);
  fclose(R2_Fastq_Ptr);
  
  return 0;
}

/// Now that we've reached coordinate on chromosome reference, look back and
/// clean up any previous positions from being tracked.
void Duplex_Consensus::cleanupPriorReads(bam1_t *b)
{
  static DupKey emptyKey;
  static DupKey tempKey2;
  
  /// Set where to stop cleaning out the structures.
  /// Initialize to the end of the structures.
  PairedMap::iterator pairedFinish = myPairedMap.end();
  //uint64_t mateStopPos = 0;
  
  /// If a record was specified, stop before this record.
  if (b != NULL)
  {
    int32_t reference = b->core.tid;
    int32_t coordinate = b->core.pos;
    tempKey2.cleanupKey(reference, coordinate);
    /// Now do the same thing with the paired reads
    PairedKey pairedKey(emptyKey, tempKey2);
    pairedFinish = myPairedMap.lower_bound(pairedKey);
    //mateStopPos = combineChromPos(reference, coordinate);
  }
  
  ///  Now do the same thing with the paired reads
  for (PairedMap::iterator piter = myPairedMap.begin();
    piter != pairedFinish; piter = myPairedMap.upper_bound(piter->first))
  {
    /// found families associated with this pair key
    /// do the consensus making here
    handleFamilyPair(piter);
    /// free the memory
    //cleanFamilyPair(piter);
  }
  
  ///  Erase the entries.
  if (pairedFinish != myPairedMap.begin())
    myPairedMap.erase(myPairedMap.begin(), pairedFinish);
  
  return;
}


/// determine whether the record's position is different from the previous record
bool Duplex_Consensus::hasPositionChanged(bam1_t *b)
{
  if (lastReference != b->core.tid ||
      lastCoordinate < b->core.pos)
  {
    if (lastReference != b->core.tid)
    {
      lastReference = b->core.tid;
    }
    lastCoordinate = b->core.pos;
    return true;
  }
  return false;
}

// When a record is read, check if it is a duplicate or
// store for future checking.
void Duplex_Consensus::checkDuplex(bam1_t *b)
{
  BamAlignment record(b);
  /// Get the key for this record.
  static DupKey key;
  key.initKey(record);
  
  int32_t chromID = b->core.tid;
  int32_t mateChromID = b->core.mtid;
  /// This is a paired record, so check for its mate.
  uint64_t readPos = combineChromPos(chromID, b->core.pos);
  uint64_t matePos = combineChromPos(mateChromID, b->core.mpos);
  
  bool foundMate = false;
 // MateData mateData;
 // MateData * mateData=new MateData;
  shared_ptr<MateData> mateData(new MateData());
  /// this is mate2, check to see if the mate is prior to this record.
  if (matePos <= readPos)
  {
    /// The mate map is stored by the mate position, so look for this
    /// record's position.
    /// The mate should be in the mate map, so find it.
    std::pair<MateMap::iterator, MateMap::iterator> matches =
      myMateMap.equal_range(readPos);
    /// Loop through the elements that matched the pos looking for the mate.
    for (MateMap::iterator iter = matches.first;
      iter != matches.second; iter++)
    {
      if (strcmp((*iter).second->readName.c_str(),
                 bam_get_qname(b)) == 0)
      {
        /// Found the matched mate1
       // std::tr1::shared_ptr<MateData> mateData1((*iter).second);
       // mateData=mateData1;
        mateData=(*iter).second;
      //  mateData = &((*iter).second);
        
        /// Remove the entry from the map.
        myMateMap.erase(iter);
        foundMate = true;
        break;
      }
    }
  }
  
  /// this is mate1
  if (!foundMate)
  {
    
    if (matePos >= readPos)
    {
      std::string barcode = get_barcode(b);
      int BC_index = get_barcode_index_bc(barcode, barcodes_list_sino);
      /// this is IDT barcode
      if (barcode.size() == 5)
      {
        barcode = trim_barcode(barcode);
        BC_index = get_barcode_index_bc_exact(barcode, barcodes_list_idt);
      }
      if (BC_index == 1000)
      {
        reads_with_wrong_barcodes++;
        return;
      }
      reads_correct_barcodes++;
      /// Haven't gotten to the mate yet, so store this record in myMateMap
      /// the key is the matePos
      // MateData mymate =new MateData(b,key,BC_index);
      shared_ptr<MateData> mateData(new MateData(b,key,BC_index,true));
      
      //mateData=mateData1;
    //  myMateMap.insert(std::make_pair(matePos, MateData(b, key, BC_index, true)));
      myMateMap.insert(std::make_pair(matePos,mateData));
      return;
    }
  }
  
  /// now the mate1 is found, make a pair between mate1 and mate2
  if (foundMate)
  {
    std::string barcode = get_barcode(b);
    int BC2_index = get_barcode_index_bc(barcode, barcodes_list_sino);
    /// this is IDT barcode
    if (barcode.size() == 5)
    {
      barcode = trim_barcode(barcode);
      BC2_index = get_barcode_index_bc_exact(barcode, barcodes_list_idt);
    }
    if (BC2_index == 1000)
    {
      reads_with_wrong_barcodes++;
      /// this is not our barcode
      return;
    }
    reads_correct_barcodes++;
    /// here the mate is found, mateData is mate1, current record b is mate2
    /// Make the paired key.
    shared_ptr<MateData> mateData1(new MateData(b,key,BC2_index,true));
    
    PairedKey pkey(mateData->key, key);
    
    PairedData pdata(mateData, mateData1);
    //std::tr1::shared_ptr<PairedData> pdata(new PairedData(mateData,mateData1));
    myPairedMap.insert(std::make_pair(pkey, pdata));
    
    read_pairs_sscs++;
    return;
  }
  else
  {
    /// here the mate1 has the wrong barcode
    n_pass_mate++;
  }
  
  return;
  
}

/**
 * @brief handle family to make consensus read pairs
 * @param familyIterator
 */
void Duplex_Consensus::handleFamilyPair(PairedMap::iterator familyIterator)
{
  std::pair<PairedMap::iterator, PairedMap::iterator> res = myPairedMap.equal_range(familyIterator->first);
  BarcodePairedMap myBarcodePairedMap;
  std::vector<PairInfo> SSfamily_vec;
  for (PairedMap::iterator i = res.first; i != res.second; ++i)
  {
    int bc1 = i->second.m1->BC_index;
    int bc2 = i->second.m2->BC_index;
    BarcodeKey barKey(bc1, bc2, i->second.m1->flag);
    PairedData pd = i->second;
    myBarcodePairedMap.insert(std::make_pair(barKey, pd));
  }
  
  for (BarcodePairedMap::iterator iter = myBarcodePairedMap.begin();
    iter != myBarcodePairedMap.end(); iter = myBarcodePairedMap.upper_bound(iter->first))
  {
    total_families++;
    /// found families associated with this pair key
    int sizeoffamily = myBarcodePairedMap.count(iter->first);
    
    map<int, int>::iterator it_s = family2size.find(sizeoffamily);
    if (it_s != family2size.end())
    {
      /// increase the counting of this family size
      (it_s->second)++;
    }
    else
    {
      family2size[sizeoffamily] = 1;
    }
    
    std::pair<BarcodePairedMap::iterator, BarcodePairedMap::iterator> res1 = myBarcodePairedMap.equal_range(
      iter->first);
    
    map<std::string, int> cigar2size;
    
    for (BarcodePairedMap::iterator i = res1.first; i != res1.second; ++i)
    {
      std::string cigar_string = i->second.m1->cigar_str + "#" + i->second.m2->cigar_str;
      map<std::string, int>::iterator it_cigar = cigar2size.find(cigar_string);
      if (it_cigar == cigar2size.end())
        cigar2size[cigar_string] = 1;
      else
        (it_cigar->second)++;
    }
    ///now find the most abundant cigar
    int cigar_count = 0;
    std::string abundant_cigar;
    
    for (map<std::string, int>::iterator it_cigar = cigar2size.begin(); it_cigar != cigar2size.end(); it_cigar++)
    {
      if (it_cigar->second > cigar_count)
      {
        abundant_cigar = std::string(it_cigar->first);
        cigar_count = it_cigar->second;
      }
    }
    
    std::vector<uint8_t *> read_seq_m1_t(cigar_count);
    std::vector<uint8_t *> read_seq_m2_t(cigar_count);
    std::vector<uint8_t *> read_qual_m1_t(cigar_count);
    std::vector<uint8_t *> read_qual_m2_t(cigar_count);
    int numreads = 0;
    
    for (BarcodePairedMap::iterator i = res1.first; i != res1.second; ++i)
    {
      std::string cigar_string = i->second.m1->cigar_str + "#" + i->second.m2->cigar_str;
      if (cigar_string == abundant_cigar)
      {
        read_seq_m1_t[numreads] = i->second.m1->seq;
        read_seq_m2_t[numreads] = i->second.m2->seq;
        read_qual_m1_t[numreads] = i->second.m1->qual;
        read_qual_m2_t[numreads] = i->second.m2->qual;
        numreads++;
      }
    }
    
    PairInfo sscs;
    int read_len = iter->second.m1->l_qseq;
    std::stringstream line;
    line << number2chrom(familyIterator->first.key1.reference) << ":" << familyIterator->first.key1.coordinate + 1
         << "#" << number2chrom(familyIterator->first.key2.reference) << ":"
         << familyIterator->first.key2.coordinate + 1
         << "#" << iter->first.bc1 + 1 << "#" << iter->first.bc2 + 1
         << "#SSCS_" << cigar_count;
    sscs.m_qname = line.str();
    sscs.m1_size = cigar_count;
    
    form_ss_consensus_phred(read_seq_m1_t, read_seq_m2_t, read_qual_m1_t, read_qual_m2_t, sscs, read_len);
    
    sscs.m1_flag = iter->second.m1->flag;
    sscs.m2_flag = iter->second.m2->flag;
    sscs.bc1 = iter->first.bc1;
    sscs.bc2 = iter->first.bc2;
    
    SSfamily_vec.push_back(sscs);
    cigar2size.erase(cigar2size.begin(), cigar2size.end());
    
    read_seq_m1_t.clear();
    read_seq_m2_t.clear();
    read_qual_m1_t.clear();
    read_qual_m2_t.clear();
    
  } //end for
  
  
  set<int> duplex_index;
  duplex_index.clear();
  set<int>::iterator set_it;
  /// go through the sscs families to find duplex
  for (int a = 0; a < SSfamily_vec.size(); a++)
  {
    sscs_families++;

    /// only do single strand consensus
    if (onlySingleStrand)
    {
      sscs_reads++;
      bool to_write = false;
    
      if (SSfamily_vec[a].m1_size > familysize_cutoff)
        to_write = true;
      else
      {
        if (considerOne)
        {
          n_ignored_family++;
          to_write = true;
        }
      }
      if (to_write)
      {
        reads_2_fastq++;
        /// can not form DS
        if (SSfamily_vec[a].m1_flag & BAM_FREAD1)
        {
          /// mate1  to R1 fastq file
          fastq_out(R1_Fastq_Ptr, SSfamily_vec[a], 1, false);
          fastq_out(R2_Fastq_Ptr, SSfamily_vec[a], 2, false);
        }
        else
        {
          /// mate1  to R2 fastq file
          fastq_out(R2_Fastq_Ptr, SSfamily_vec[a], 1, false);
          fastq_out(R1_Fastq_Ptr, SSfamily_vec[a], 2, false);
        }
      
      }
      else
        n_ignored_family++;
    
    }
  
  
    //// do double strand consensus
    else
    {
      bool foundDS = false;
      int DS_index = -1;
      
      for (int b = a + 1; b < SSfamily_vec.size(); b++)
      {
        set_it = duplex_index.find(b);
        if (set_it != duplex_index.end())
          continue;
      
        if (isDuplex(SSfamily_vec[a], SSfamily_vec[b]))
        {
          foundDS = true;
          DS_index = b;
          break;
        }
      }
      /// form the duplex
    
      if (foundDS)
      {
        dscs_reads++;
        duplex_index.insert(DS_index);
        duplex_index.insert(a);
        PairInfo dscs;
        dscs.m1_size = SSfamily_vec[a].m1_size;
        dscs.m2_size = SSfamily_vec[DS_index].m1_size;
        bool ds_fq;
        /// for those dscs with one top and one bottom family
        //if (dscs.m1_size == 1 && dscs.m1_size == dscs.m2_size)
        // ds_fq = generate_ds_consensus_seq_phred(SSfamily_vec[a], SSfamily_vec[DS_index], dscs);
        //else
        ds_fq = generate_ds_consensus_seq(SSfamily_vec[a], SSfamily_vec[DS_index], dscs);
      
        if (ds_fq)
        {
          reads_2_fastq++;
          if (dscs.m1_flag & BAM_FREAD1)
          {
            /// mate1  to R1 fastq file
            fastq_out(R1_Fastq_Ptr, dscs, 1, false);
            fastq_out(R2_Fastq_Ptr, dscs, 2, false);
          }
          else
          {
            /// mate1  to R2 fastq file
            fastq_out(R2_Fastq_Ptr, dscs, 1, false);
            fastq_out(R1_Fastq_Ptr, dscs, 2, false);
          }
        }
      
      }
      else
      {
        set_it = duplex_index.find(a);
        if (set_it == duplex_index.end())
        {
          sscs_reads++;
          bool to_write = false;
        
          if (SSfamily_vec[a].m1_size > familysize_cutoff)
            to_write = true;
          else
          {
            if (considerOne)
            {
              n_ignored_family++;
              to_write = true;
            }
          }
          if (to_write)
          {
            reads_2_fastq++;
            /// can not form DS
            if (SSfamily_vec[a].m1_flag & BAM_FREAD1)
            {
              /// mate1  to R1 fastq file
              fastq_out(R1_Fastq_Ptr, SSfamily_vec[a], 1, false);
              fastq_out(R2_Fastq_Ptr, SSfamily_vec[a], 2, false);
            }
            else
            {
              /// mate1  to R2 fastq file
              fastq_out(R2_Fastq_Ptr, SSfamily_vec[a], 1, false);
              fastq_out(R1_Fastq_Ptr, SSfamily_vec[a], 2, false);
            }
          
          }
          else
            n_ignored_family++;
        }
      }
    }
  
  }
  
  duplex_index.clear();
  SSfamily_vec.erase(SSfamily_vec.begin(), SSfamily_vec.end());
  myBarcodePairedMap.erase(myBarcodePairedMap.begin(), myBarcodePairedMap.end());
}


bool Duplex_Consensus::form_ss_consensus(std::vector<uint8_t *> read_seq_m1,
                                         std::vector<uint8_t *> read_seq_m2,
                                         std::vector<uint8_t *> read_qual_m1,
                                         std::vector<uint8_t *> read_qual_m2,
                                         PairInfo &sscs, int read_length)
{
  //int family_size_cigar = read_seq_m1.size();
  int m1_N_count = 0;
  int m2_N_count = 0;
  
  for (int j = 0; j < read_length; j++)
  {
    int countNucl[5][2] = {{0}};
    int m1_count = 0;
    int m2_count = 0;
    
    for (size_t jj = 0; jj < read_seq_m1.size(); jj++)
    {
      
      if (read_qual_m1[jj][j] >= myMinQual)
      {
        int base_m1 = bam_seqi(read_seq_m1[jj], j);
        /// here is base "A"
        if (base_m1 == 1)
          countNucl[0][0] = countNucl[0][0] + 1;
        /// base "C"
        if (base_m1 == 2)
          countNucl[1][0] = countNucl[1][0] + 1;
        /// base "G"
        if (base_m1 == 4)
          countNucl[2][0] = countNucl[2][0] + 1;
        /// base "T"
        if (base_m1 == 8)
          countNucl[3][0] = countNucl[3][0] + 1;
        /// base "N"
        if (base_m1 == 15)
          countNucl[4][0] = countNucl[4][0] + 1;
        m1_count++;
      }
      
      //mate2
      if (read_qual_m2[jj][j] >= myMinQual)
      {
        int base_m2 = bam_seqi(read_seq_m2[jj], j);
        if (base_m2 == 1)
          countNucl[0][1] = countNucl[0][1] + 1;
        if (base_m2 == 2)
          countNucl[1][1] = countNucl[1][1] + 1;
        if (base_m2 == 4)
          countNucl[2][1] = countNucl[2][1] + 1;
        if (base_m2 == 8)
          countNucl[3][1] = countNucl[3][1] + 1;
        if (base_m2 == 15)
          countNucl[4][1] = countNucl[4][1] + 1;
        m2_count++;
      }
    }
    
    int maxNucl_m1 = countNucl[0][0];
    string consNucl_m1 = Nucl[0];
    
    int maxNucl_m2 = countNucl[0][1];
    string consNucl_m2 = Nucl[0];
    
    string frscore_m1;
    string frscore_m2;
    
    for (int kk = 1; kk < 5; kk++)
    {
      if (countNucl[kk][0] > maxNucl_m1)
      {
        maxNucl_m1 = countNucl[kk][0];
        consNucl_m1 = Nucl[kk];
      }
      
      if (countNucl[kk][1] > maxNucl_m2)
      {
        maxNucl_m2 = countNucl[kk][1];
        consNucl_m2 = Nucl[kk];
      }
    }
    
    if ((maxNucl_m1 / (double) m1_count) >= abundant_nucleotide_cutoff)
    {
      sscs.m1_seq.append(consNucl_m1);
      frscore_m1 = (char) (frscore_sscs + phred_score);
      sscs.m1_qual.append(frscore_m1);
    }
    else
    {
      m1_N_count++;
      sscs.m1_seq.append("N");
      frscore_m1 = (char) (frscore_N + phred_score);
      sscs.m1_qual.append(frscore_m1);
    }
    
    if ((maxNucl_m2 / (double) m2_count) >= abundant_nucleotide_cutoff)
    {
      sscs.m2_seq.append(consNucl_m2);
      frscore_m2 = (char) (frscore_sscs + phred_score);
      sscs.m2_qual.append(frscore_m2);
    }
    else
    {
      m2_N_count++;
      sscs.m2_seq.append("N");
      frscore_m2 = (char) (frscore_N + phred_score);
      sscs.m2_qual.append(frscore_m2);
    }
    
  }
  
  float N_ratio1 = (float) (m1_N_count) / (float) (read_length);
  float N_ratio2 = (float) (m2_N_count) / (float) (read_length);
  map<int, int>::iterator it_s = DS_N_size.find(m1_N_count);
  if (it_s != DS_N_size.end())
    it_s->second++;
  else
    DS_N_size[m1_N_count] = 1;
  
  it_s = DS_N_size.find(m2_N_count);
  if (it_s != DS_N_size.end())
    it_s->second++;
  else
    DS_N_size[m2_N_count] = 1;
  
  if ((N_ratio1 >= N_percent_cutoff) || (N_ratio2 >= N_percent_cutoff))
  {
    N_SS_family++;
    return false;
  }
  else
    return true;
}


void Duplex_Consensus::form_ss_consensus_phred(std::vector<uint8_t *> read_seq_m1,
                                               std::vector<uint8_t *> read_seq_m2,
                                               std::vector<uint8_t *> read_qual_m1,
                                               std::vector<uint8_t *> read_qual_m2,
                                               PairInfo &sscs, int read_length)
{
  int family_size = read_seq_m1.size();
  for (int j = 0; j < read_length; j++)
  {
    int m1_count = 0;
    int m2_count = 0;
    double_t error4_m1[4][family_size];
    double_t error4_m2[4][family_size];
    
    for (int jj = 0; jj < family_size; jj++)
    {
      int base_m1 = bam_seqi(read_seq_m1[jj], j);
      int baseQ_m1 = read_qual_m1[jj][j];
      
      for (int i = 0; i < 4; ++i)
      {
        error4_m1[i][m1_count] = get_error_pr(base5[i], base_m1, baseQ_m1);
      }
      m1_count++;
      
      
      int base_m2 = bam_seqi(read_seq_m2[jj], j);
      int baseQ_m2 = read_qual_m2[jj][j];
      
      for (int i = 0; i < 4; ++i)
      {
        error4_m2[i][m2_count] = get_error_pr(base5[i], base_m2, baseQ_m2);
      }
      m2_count++;
      
    }
    
    double_t error4_m1_sum[4] = {0.0};
    double_t error4_m2_sum[4] = {0.0};
    /// sum up the error probabilities
    for (int i = 0; i < 4; ++i)
    {
      error4_m1_sum[i] = error4_m1[i][0];
      error4_m2_sum[i] = error4_m2[i][0];
      for (int k = 1; k < m1_count; k++)
      {
        error4_m1_sum[i] *= error4_m1[i][k];
      }
      for (int k = 1; k < m2_count; k++)
      {
        error4_m2_sum[i] *= error4_m2[i][k];
      }
    }
    
    int max_base_m1 = 0;
    double_t max_score_m1 = error4_m1_sum[0];
    int max_base_m2 = 0;
    double_t max_score_m2 = error4_m2_sum[0];
    
    double_t sum_m1 = 0.0;
    double_t sum_m2 = 0.0;
    for (int i = 0; i < 4; i++)
    {
      sum_m1 += error4_m1_sum[i];
      sum_m2 += error4_m2_sum[i];
      if (error4_m1_sum[i] > max_score_m1)
      {
        max_base_m1 = i;
        max_score_m1 = error4_m1_sum[i];
      }
      if (error4_m2_sum[i] > max_score_m2)
      {
        max_base_m2 = i;
        max_score_m2 = error4_m2_sum[i];
      }
    }
    int base_phred_m1, base_phred_m2;
    if (max_score_m1 == sum_m1)
      base_phred_m1 = 93;
    else
    {
      base_phred_m1 = int(-10 * log10(1 - max_score_m1 / sum_m1));
      if (base_phred_m1 >= 93)
        base_phred_m1 = 93;
    }
    if (max_score_m2 == sum_m2)
      base_phred_m2 = 93;
    else
    {
      base_phred_m2 = int(-10 * log10(1 - max_score_m2 / sum_m2));
      if (base_phred_m2 >= 93)
        base_phred_m2 = 93;
    }
    
    
    /// assign the base and its quality
    sscs.m1_seq += bam_nt16_rev_table[base5[max_base_m1]];
    sscs.m1_qual += (char) (base_phred_m1 + phred_score);
    
    
    sscs.m2_seq += bam_nt16_rev_table[base5[max_base_m2]];
    sscs.m2_qual += (char) (base_phred_m2 + phred_score);
    
  }
  
  return;
}

void Duplex_Consensus::fastq_out(FILE *f, PairInfo &read_pair, int mate, bool append12)
{
  
  kstring_t linebuf = {0, 0, NULL};
  kputsn("", 0, &linebuf);
  
  linebuf.l = 0;
  /// write read name
  kputc('@', &linebuf);
  kputs((read_pair.m_qname).c_str(), &linebuf);
  
  // mate 1
  if (mate == 1)
  {
    /// write sequence
    if (append12) kputs("/1\n", &linebuf);
    else kputc('\n', &linebuf);
    
    if (read_pair.m1_flag & BAM_FREVERSE)
    {
      reverseComplement(read_pair.m1_seq);
      reverseSequence(read_pair.m1_qual);
    }
    
    kputs(read_pair.m1_seq.c_str(), &linebuf);
    kputc('\n', &linebuf);
    
    kputs("+\n", &linebuf);
    /// write quality string
    kputs(read_pair.m1_qual.c_str(), &linebuf);
    kputc('\n', &linebuf);
  }
    /// mate 2
  else
  {
    if (append12) kputs("/2\n", &linebuf);
    else kputc('\n', &linebuf);
    /// write sequence
    
    if (read_pair.m2_flag & BAM_FREVERSE)
    {
      reverseComplement(read_pair.m2_seq);
      reverseSequence(read_pair.m2_qual);
    }
    
    kputs(read_pair.m2_seq.c_str(), &linebuf);
    kputc('\n', &linebuf);
    
    kputs("+\n", &linebuf);
    /// write quality
    kputs(read_pair.m2_qual.c_str(), &linebuf);
    kputc('\n', &linebuf);
    
  }
  fputs(linebuf.s, f);
  free(linebuf.s);
}

/**
 * @brief Generate double strand consensus sequence and quality strings from two single strand consensus
 *
 * @param top top pair of a SSCS
 * @param bottom bottom pair of a SSCS
 * @param[out] dscs output double strand consensus sequence
 */
bool Duplex_Consensus::generate_ds_consensus_seq(PairInfo &top, PairInfo &bottom, PairInfo &dscs)
{
  string frscore_m1, frscore_m2;
  string m1_nucl, m2_nucl;
  int m1_N_count = 0;
  int m2_N_count = 0;
  dscs.m1_flag = top.m1_flag;
  dscs.m2_flag = top.m2_flag;
  
  
  vector<string> b_ss = split_string(bottom.m_qname, "#");
  string bottom_info = b_ss[b_ss.size() - 1];
  
  dscs.m_qname = top.m_qname + "#" + bottom_info;
  
  string m1_seq_top = top.m1_seq;
  string m2_seq_top = top.m2_seq;
  
  string m1_seq_bottom = bottom.m1_seq;
  string m2_seq_bottom = bottom.m2_seq;
  
  int read_length = m1_seq_top.size();
  
  for (int j = 0; j < read_length; j++)
  {
    //Create DSCS for mate1
    
    //If the nucleotide is the same in both sscs
    //then the nucleotide is written into the dssc.
    //Otherwise that position is marked with an N
    if (m1_seq_top[j] == m1_seq_bottom[j])
    {
      m1_nucl = m1_seq_top[j];
      dscs.m1_seq.append(m1_nucl); //nucleotide
      frscore_m1 = (char) (frscore_dscs + phred_score);
      dscs.m1_qual.append(frscore_m1);//phred_score of the nucleotide
  
    }
    else
    {
      dscs.m1_seq.append("N");
      frscore_m1 = (char) (frscore_N_dscs + phred_score);
      dscs.m1_qual.append(frscore_m1);
      m1_N_count++;
    }
    
    /// Create DSCS for mate2
    if (m2_seq_top[j] == m2_seq_bottom[j])
    {
      m2_nucl = m2_seq_top[j];
      dscs.m2_seq.append(m2_nucl);
      frscore_m2 = (char) (frscore_dscs + phred_score);
      dscs.m2_qual.append(frscore_m2);
    }
    else
    {
      dscs.m2_seq.append("N");
      frscore_m2 = (char) (frscore_N_dscs + phred_score);
      dscs.m2_qual.append(frscore_m2);
      m2_N_count++;
    }
    
  }
  float N_ratio1 = (float) (m1_N_count) / (float) (read_length);
  float N_ratio2 = (float) (m2_N_count) / (float) (read_length);
  map<int, int>::iterator it_s = DS_N_size.find(m1_N_count);
  if (it_s != DS_N_size.end())
    it_s->second++;
  else
    DS_N_size[m1_N_count] = 1;
  
  it_s = DS_N_size.find(m2_N_count);
  if (it_s != DS_N_size.end())
    it_s->second++;
  else
    DS_N_size[m2_N_count] = 1;
  
  if ((N_ratio1 >= N_percent_cutoff) || (N_ratio2 >= N_percent_cutoff))
  {
    N_DS_family++;
//    cout<<top.m1_seq<<endl;
//    cout<<bottom.m1_seq<<endl;
//    cout<<dscs.m1_seq<<endl;
//    cout<<endl;
    return false;
  }
  else
    return true;
}

bool Duplex_Consensus::generate_ds_consensus_seq_phred(PairInfo &top, PairInfo &bottom, PairInfo &dscs)
{
  string frscore_m1, frscore_m2;
  string m1_nucl, m2_nucl;
  int m1_N_count = 0;
  int m2_N_count = 0;
  dscs.m1_flag = top.m1_flag;
  dscs.m2_flag = top.m2_flag;
  vector<string> b_ss = split_string(bottom.m_qname, "#");
  string bottom_info = b_ss[b_ss.size() - 1];
  
  dscs.m_qname = top.m_qname + "#" + bottom_info;
  
  string m1_seq_top = top.m1_seq;
  string m2_seq_top = top.m2_seq;
  
  string m1_seq_bottom = bottom.m1_seq;
  string m2_seq_bottom = bottom.m2_seq;
  
  int read_length = m1_seq_top.size();
  
  for (int j = 0; j < read_length; j++)
  {
    //Create DSCS for mate1
    
    //If the nucleotide is the same in both sscs
    //then the nucleotide is written into the dssc.
    //Otherwise that position is marked with an N
    if (m1_seq_top[j] == m1_seq_bottom[j])
    {
      m1_nucl = m1_seq_top[j];
      dscs.m1_seq.append(m1_nucl); //nucleotide
      int bq = bottom.m1_qual[j] - phred_score;
      int tq = top.m1_qual[j] - phred_score;
      double bq_m1 = 1 - pow(10, -0.1 * bq);
      double tq_m1 = 1 - pow(10, -0.1 * tq);
      double score_m1 = bq_m1 * tq_m1;
      double sum_er = bq_m1 * tq_m1 + pow(10, -0.1 * bq) * pow(10, -0.1 * bq)
                      + pow(10, -0.1 * tq) * pow(10, -0.1 * tq);
      
      int base_phred_m1;
      if (score_m1 == sum_er)
        base_phred_m1 = 93;
      else
      {
        base_phred_m1 = int(-10 * log10(1 - score_m1 / sum_er));
        if (base_phred_m1 > 93)
          base_phred_m1 = 93;
      }
      
      frscore_m1 = (char) (base_phred_m1 + phred_score);
      dscs.m1_qual.append(frscore_m1);//phred_score of the nucleotide
      
    }
    else
    {
      dscs.m1_seq.append("N");
      frscore_m1 = (char) (frscore_N_dscs + phred_score);
      dscs.m1_qual.append(frscore_m1);
      m1_N_count++;
    }
    
    /// Create DSCS for mate2
    if (m2_seq_top[j] == m2_seq_bottom[j])
    {
      m2_nucl = m2_seq_top[j];
      dscs.m2_seq.append(m2_nucl);
      int bq = bottom.m2_qual[j] - phred_score;
      int tq = top.m2_qual[j] - phred_score;
      double bq_m2 = 1 - pow(10, -0.1 * bq);
      double tq_m2 = 1 - pow(10, -0.1 * tq);
      double score_m2 = bq_m2 * tq_m2;
      double sum_er = bq_m2 * tq_m2 + pow(10, -0.1 * bq) * pow(10, -0.1 * bq)
                      + pow(10, -0.1 * tq) * pow(10, -0.1 * tq);
      
      int base_phred_m2;
      if (score_m2 == sum_er)
        base_phred_m2 = 93;
      else
      {
        base_phred_m2 = int(-10 * log10(1 - score_m2 / sum_er));
        if (base_phred_m2 > 93)
          base_phred_m2 = 93;
      }
      frscore_m2 = (char) (base_phred_m2 + phred_score);
      dscs.m2_qual.append(frscore_m2);
      
    }
    else
    {
      dscs.m2_seq.append("N");
      frscore_m2 = (char) (frscore_N_dscs + phred_score);
      dscs.m2_qual.append(frscore_m2);
      m2_N_count++;
    }
    
  }
  float N_ratio1 = (float) (m1_N_count) / (float) (read_length);
  float N_ratio2 = (float) (m2_N_count) / (float) (read_length);
  map<int, int>::iterator it_s = DS_N_size.find(m1_N_count);
  if (it_s != DS_N_size.end())
    it_s->second++;
  else
    DS_N_size[m1_N_count] = 1;
  
  it_s = DS_N_size.find(m2_N_count);
  if (it_s != DS_N_size.end())
    it_s->second++;
  else
    DS_N_size[m2_N_count] = 1;
  
  if ((N_ratio1 >= N_percent_cutoff) || (N_ratio2 >= N_percent_cutoff))
  {
    N_DS_family++;
//    cout<<top.m1_seq<<endl;
//    cout<<bottom.m1_seq<<endl;
//    cout<<dscs.m1_seq<<endl;
//    cout<<endl;
    return false;
  }
  else
    return true;
}

/**
 * @brief use flat to get check these two pairs for a duplex
 * @param p1
 * @param p2
 * @return
 */
bool Duplex_Consensus::isDuplex(PairInfo &p1, PairInfo &p2)
{
  bool duplex = false;
  bool p1_R1 = p1.m1_flag & BAM_FREAD1;
  bool p2_R1 = p2.m1_flag & BAM_FREAD1;
  /// if P1.m1 is R1 and p2.m1 is R2  or
  /// if P1.m1 is R2 and p2.m1 is R1
  if ((p1_R1 && !(p2_R1)) || (!(p1_R1) && p2_R1))
    duplex = true;
  
  return duplex;
}

/**
 * @brief clean up the map to free memory
 * @param familyIterator
 */
void Duplex_Consensus::cleanFamilyPair(PairedMap::iterator familyIterator)
{
  std::pair<PairedMap::iterator, PairedMap::iterator> res = myPairedMap.equal_range(familyIterator->first);
  for (PairedMap::iterator i = res.first; i != res.second; ++i)
  {
    if (i->second.m1->seq != NULL)
    {
      free(i->second.m1->seq);
      i->second.m1->seq = NULL;
    }
    if (i->second.m1->qual != NULL)
    {
      free(i->second.m1->qual);
      i->second.m1->qual = NULL;
    }
    
    if (i->second.m2->seq != NULL)
    {
      free(i->second.m2->seq);
      i->second.m2->seq = NULL;
    }
    if (i->second.m2->qual != NULL)
    {
      free(i->second.m2->qual);
      i->second.m2->qual = NULL;
    }
  }
  
}

double_t Duplex_Consensus::get_error_pr(int givenBase, int base, int baseQ)
{
  double_t error_pr;
  if (base == givenBase)
  {
    error_pr = 1 - pow(10, -0.1 * (double) baseQ);
  }
  else
  {
    error_pr = pow(10, -0.1 * (double) baseQ) / 3.0;
  }
  
  return error_pr;
  
}