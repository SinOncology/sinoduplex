#pragma once
#include <iostream>

// Logger class for logging/error/warning
class Logger {
 protected:
  FILE* fp_log;
  FILE* fp_err;
  bool b_verbose;

  Logger() {} // default constructor prohibited
 public:
  static Logger* gLogger;
  Logger(const char* filename, bool verbose);
  void writeLog(const char* format, ...);
  void error(const char* format, ...);
  void warning(const char* format, ...);
};

