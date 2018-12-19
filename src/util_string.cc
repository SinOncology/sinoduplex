
#include "util_string.h"


/**
 * @brief Performs an in-place sequence reversal
 * 
 * @param[in,out] seq the sequence to reverse
 */
void reverseSequence(std::string &seq)
{
  reverse(seq.begin(), seq.end());
}

/**
 *@brief Performs an in-place reverse complement conversion
 * 
 * @param[in,out] seq the sequence to reverse complement
 */
void reverseComplement(std::string &seq)
{
  
  // reverse the sequence
  reverseSequence(seq);
  
  // swap the bases
  for (unsigned int i = 0; i < seq.length(); i++)
  {
    switch (seq[i])
    {
    case 'A':
      seq[i] = 'T';
      break;
    case 'C':
      seq[i] = 'G';
      break;
    case 'G':
      seq[i] = 'C';
      break;
    case 'T':
      seq[i] = 'A';
      break;
    case 'a':
      seq[i] = 't';
      break;
    case 'c':
      seq[i] = 'g';
      break;
    case 'g':
      seq[i] = 'c';
      break;
    case 't':
      seq[i] = 'a';
      break;
    default:
      break;
    }
  }
}


std::vector <std::string> split_string(const std::string &s,
                                       const std::string &delim)
{
  bool keep_empty = false;
  std::vector <std::string> result;
  result.empty();
  if (delim.empty())
  {
    result.push_back(s);
    return result;
  }
  
  std::string::const_iterator substart = s.begin(), subend;
  while (true)
  {
    subend = search(substart, s.end(), delim.begin(), delim.end());
    std::string temp(substart, subend);
    if (keep_empty || !temp.empty())
    {
      result.push_back(temp);
    }
    if (subend == s.end())
    {
      break;
    }
    substart = subend + delim.size();
  }
  return result;
}




/*
 * usage: judge the input string is fully composed by repeat unit
 */
pair<int, string> strRepeatEasy(const string &str)
{
  //violent window
  int repeatCount;
  // bool repeat = false;
  std::string sub = "", current_sub;
  unsigned long len = str.size();
  pair<int, string> result_pair;
  for (int i = 1; i < str.size() / 2; ++i)
  {
    //i : window size
    if (len % i != 0) continue;//if window size is not exact division of string size : illegal
    sub = str.substr(0, i);
    repeatCount = 1;
    for (int j = 1; j <= str.size() / i; ++j)
    {
      //j: repeat unit count
      current_sub = str.substr(j * i, i);
      if (current_sub != sub)
      {
        break;
      }
      else
      {
        repeatCount++;
      }
    }
    if (repeatCount == str.size() / i)
    {
      break;
    }
    else
    {
      repeatCount = 0;
      sub = "";
    }
  }
  return make_pair(repeatCount, sub);
}

int countCharInStr(string str, char i)
{
  return 0;
}



string condense_repeat_string(string expanded_string)
{
  string condense = "";
  int count = 1;
  char c_c;
  //char c_n;
  char c_p = expanded_string[0];
  for (int i = 1; i < expanded_string.size();)
  {
    c_c = expanded_string[i];
    if (c_c == c_p)
    {
      c_p = c_c;
      count++;
      i++;
    }
    else
    {
      condense += to_string(count) + c_p;
      count = 1;
      c_p = c_c;
      i++;
    }
    
  }
  condense += to_string(count) + c_p;
  return condense;
}



/**
 * @brief check a given in the barcode list and return the index of
 * BC in the barcode list
 * @param bc
 * @param barcodes_list
 * @return
 */
int get_barcode_index_bc(string bc, const std::vector<std::string> &barcodes_list)
{
  int count_BC;
  int BC_index;
  int index_BC;
  
  /// Maximum number of mismatches allowed to accept a barcode
  /// here is 1 mismatch
  const int max_nmis = 1;
  int min_mis_BC = 1000;
  index_BC = 1000;
  
  /// Determine if the read contains a good barcode or not
  for (std::size_t j = 0; j < barcodes_list.size(); j++)
  {
    count_BC = 0;
    for (std::size_t k = 0; k < barcodes_list[j].size(); k++)
    {
      /// Count mismatches for BC mate1
      if (barcodes_list[j][k] != bc[k])
        count_BC++;
    }
    if (count_BC < min_mis_BC)
    {
      min_mis_BC = count_BC;
      index_BC = j;
    }
  }
  
  if (min_mis_BC <= max_nmis)
    BC_index = index_BC;
  else
    BC_index = 1000;
  
  return BC_index;
}

int get_barcode_index_bc_exact(string bc, const std::vector<std::string> &barcodes_list)
{
  int bc_index = 1000;
  
  /// exact match
  for (std::size_t j = 0; j < barcodes_list.size(); j++)
  {
    if (barcodes_list[j] == bc)
    {
      bc_index = j;
      break;
    }
  }
  return bc_index;
}

void read_barcode_list(string fn, std::vector<std::string> &barcodes_list)
{
  ifstream in;
  barcodes_list.clear();
  
  in.open(fn);
  if (!in.is_open())
  {
    cerr << "Error: cannot open barcode table " << fn << "\n";
    exit(1);
  }
  /// Number of barcodes read from the number of non-empty lines of the barcode file
  std::string current_line;
  while (std::getline(in, current_line))
  {
    current_line.erase(current_line.find_last_not_of(" \n\r\t") + 1);
    if (current_line.empty()) continue;
    barcodes_list.push_back(current_line);
  }
  in.close();
}

/**
 * @brief trim IDT barcode according to NNNT and NNNG/CT
 * @param barcode
 * @return
 */
string trim_barcode(string barcode)
{
  string trim_str;
  int len = barcode.size();
  /// this  is NNNG/CT barcode
  if (barcode.substr(len - 1, 1) == "T" && (barcode.substr(len - 2, 1) == "G" || barcode.substr(len - 2, 1) == "C"))
  {
    /// trim the last T
    trim_str = barcode.substr(0, len - 1);
  }
  else if (barcode.substr(len - 2, 1) == "T")
  {
    /// trim the last  base and T in the fourth position
    trim_str = barcode.substr(0, len - 2);
  }
  else
  {
    /// do not trim
    trim_str = barcode;
  }
  
  return trim_str;
  
}