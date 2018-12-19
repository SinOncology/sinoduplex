#pragma once
#include <string>
#include <iostream>


/// Base Class BAM Executable.
class SinoExecutable
{
public:
  static void printSinoToolsVersion(std::ostream &os);
  
  static void printSinoToolsDescription(std::ostream &os);
  
  SinoExecutable();
  
  virtual ~SinoExecutable();
  
  /// Print the
  virtual void printDescription(std::ostream &os);
  
  virtual void printUsage(std::ostream &os);
  
  virtual int execute(int argc, char **argv) = 0;
  
protected:

private:
};


