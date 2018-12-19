#include <stdlib.h>
#include <iostream>
#include "sinoduplex.h"
#include "trimFq.h"
#include "SinoExecutable.h"
#include "Duplex_Consensus.h"

std::string sinotools_usage = "     SYNOPSIS \n \t sinotools [command] <options> \n\n \
     availabl commands:\n\n \
     \t duplex      \t new approach to do duplex consensus \\n \
     \t trim_fq     \t trim fastq seq's IDT barcode\n";


SinoExecutable *CreateSinoToolsExe(const std::string &name)
{
  SinoExecutable *ret = NULL;
  

 if (name == "duplex")
  {
    ret = new Duplex_Consensus();
  }
  
  return ret;
  
}


int main(int argc, char *argv[])
{
  if (argc == 1)
  {
    print_header();
    std::cerr << sinotools_usage;
    exit(0);
  }
  
  std::string command = (std::string) argv[1];
  
  
  if (command == "duplex")
  {
    //std::cout << "duplex" << std::endl;
  }
  

  else if (command == "trim_fq")
    trim_fq(argc, argv);

  else
  {
    print_header();
    std::cerr << sinotools_usage;
    exit(0);
  }
  
  SinoExecutable *bamExe = CreateSinoToolsExe(command);
  if (bamExe != NULL)
  {
    bamExe->printDescription(std::cerr);
    int ret;
    ret = bamExe->execute(argc, argv);
    delete bamExe;
    bamExe = NULL;
  }
  
}


