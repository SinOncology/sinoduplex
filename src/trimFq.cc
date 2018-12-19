#include <zlib.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <getopt.h>
#include "htslib/kseq.h"
#include "htslib/kstring.h"

#include "trimFq.h"


using namespace std;
KSEQ_INIT(gzFile, gzread)

string trim_fq_help = "      SYNOPSIS \n \t sinotools trim_fq <options>\n\n \
     DESCRIPTION\n \
     \t -h -? -help \t help\n \
     \t -i \t input fastq file \n \
     \t -o \t output fastq file \n \
     \t -n  \t barcode_length [8] \n \
     \t -tail_n \t number of bases to be trimmed from tail [default no] \n";


void trim_fq(int argc, char *argv[])
{
  
  int longindex, opt;
  // option definition
  string fq = "";
  string output_fq = "";
  /// add option to trim from tail
  bool trim_tail = false;
  /// number of bases to be trimmed from tail
  int tail_n = 0;
  static struct option longopts[] = {
    {"help",   0, 0, 'h'},
    {"help",   0, 0, '?'},
    {"i",      1, 0, 1},
    {"o",      1, 0, 2},
    {"n",      1, 0, 3},
    {"tail_n", 1, 0, 4},
    {0,        0, 0, 0}
  };
  
  int barcode_num = 8;
  
  optind = 0;
  //parse command line arguments
  while ((opt = getopt_long_only(argc, argv, "h?", longopts, &longindex)) != -1)
  {
    switch (opt)
    {
    case 'h':
      print_header();
      cerr << trim_fq_help;
      exit(1);
    case '?':
      print_header();
      cerr << trim_fq_help;
      exit(1);
    case 1:
      fq = (string) (optarg);
      break;
    case 2:
      output_fq = (string) (optarg);
      break;
    case 3:
      barcode_num = (int) abs(atol(optarg));
      break;
    case 4:
      tail_n = (int) abs(atol(optarg));
      trim_tail = true;
      break;
    default:
      cerr << "Error: cannot parse arguments.\n";
      exit(1);
    }
  }
  if ((argc - optind) != 1)
  {
    print_header();
    cerr << trim_fq_help;
    exit(1);
  }
  
  FILE *fastq_out;
  fastq_out = fopen(output_fq.c_str(), "w");
  
  gzFile fp = NULL;
  kseq_t *seq;
  int l;
  const char *mode = "r";
  //fp = gzdopen(fileno(stdin), mode);
  fp = gzopen(fq.c_str(), mode);
  if (fp == NULL)
  {
    cerr << "can not open fastq file\n"
         << fq << endl;
    exit(1);
  }
  cout << "trim input fastq file by " << barcode_num << " bases" << endl;
  cout << "input fastq file\t" << fq << endl;
  cout << "output fastq file\t" << output_fq << endl;
  
  seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0)
  {
  
  
    kstring_t linebuf = {0, 0, NULL};
    kputsn("", 0, &linebuf);
    linebuf.l = 0;
    // write read name
    kputc('@', &linebuf);
    kputs(seq->name.s, &linebuf);
    kputs(" BC:Z:", &linebuf);
  
    kputsn(seq->seq.s, barcode_num, &linebuf);
  
    kputc('\n', &linebuf);
    if (trim_tail)
    {
      kputsn(seq->seq.s + barcode_num, (strlen(seq->seq.s) - barcode_num - tail_n), &linebuf);
    }
    else
    {
      kputs(seq->seq.s + barcode_num, &linebuf);
    }
    kputc('\n', &linebuf);
    kputs("+\n", &linebuf);
    if (trim_tail)
    {
      kputsn(seq->qual.s + barcode_num, (strlen(seq->qual.s) - barcode_num - tail_n), &linebuf);
    }
    else
    {
      kputs(seq->qual.s + barcode_num, &linebuf);
    }
    kputc('\n', &linebuf);
  
    fputs(linebuf.s, fastq_out);
  
    free(linebuf.s);
  }
  gzclose(fp);
  kseq_destroy(seq);
  fclose(fastq_out);
}
