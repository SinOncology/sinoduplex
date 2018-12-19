/**
 * @file BamExecutable.cc
 * @brief A bass class for an Executable module
 * @author Bingding Huang
 * @date 2017-01-02
 */


#include "SinoExecutable.h"


using namespace std;

SinoExecutable::SinoExecutable()
{
}


SinoExecutable::~SinoExecutable()
{
}


void SinoExecutable::printSinoToolsVersion(std::ostream &os)
{
  
  os << "Version: " <<"1.0.0"<< std::endl;
  
}

void SinoExecutable::printSinoToolsDescription(std::ostream &os)
{
  os << " --- CaGe-A: Cancer Genome Scanner Analysis --- " << std::endl;
  printSinoToolsVersion(os);
}


void SinoExecutable::printDescription(std::ostream &os)
{
  printSinoToolsDescription(os);
}


void SinoExecutable::printUsage(std::ostream &os)
{
  printSinoToolsVersion(os);
  os << endl;
  printDescription(os);
}
