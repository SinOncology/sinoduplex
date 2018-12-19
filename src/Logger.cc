
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdexcept>
#include "Logger.h"

Logger* Logger::gLogger = NULL;

// Constructor of logger
Logger::Logger(const char* filename, bool verbose) 
{
  b_verbose = verbose;
  if((filename != NULL) && (filename[0] == '-'))
  {
      // write to stderr.
      fp_log = stderr;
  }
  else
  {
      fp_log = fopen(filename, "w");
  }
  if ( fp_log == NULL ) {
      std::string errorMsg = "ERROR: Cannot open the log file ";
      errorMsg += filename;
      errorMsg += ". Check if the directory exists and you have the permission to create a file";
      fprintf(stderr,"%s",errorMsg.c_str());
      throw(std::runtime_error(errorMsg.c_str()));
  }
  fp_err = stderr;
}

// Write a log to output file
void Logger::writeLog(const char* format, ... ) {
  va_list args;
  va_start (args, format);
  vfprintf(fp_log, format, args);
  va_end (args);
  fprintf(fp_log, "\n");
  fflush(fp_log);
  
  if ( b_verbose ) {
    va_start (args, format);
    vfprintf(fp_err, format, args);
    va_end (args);
    fprintf(fp_err, "\n");
  }
}

// Write error messages and throw an exception.
void Logger::error(const char* format, ... ) {
  va_list args;
  va_start (args, format);
  fprintf(fp_log, "ERROR: ");
  vfprintf(fp_log, format, args);
  va_end (args);
  fprintf(fp_log, "\n");
  fflush(fp_log);

  va_start (args, format);
  fprintf(fp_err, "ERROR : ");
  vfprintf(fp_err, format, args);
  va_end (args);
  fprintf(fp_err, "\n");

  char buffer[256];
  va_start (args, format);
  vsnprintf (buffer,256,format, args);
  va_end (args);

  std::string errorMsg = "ERROR: ";
  errorMsg += buffer;
  throw(std::runtime_error(errorMsg.c_str()));
}

// Write warning messages
void Logger::warning(const char* format, ... ) {
  va_list args;
  va_start (args, format);
  fprintf(fp_log, "WARNING: ");
  vfprintf(fp_log, format, args);
  va_end (args);
  fprintf(fp_log, "\n");
  fflush(fp_log);

  if ( b_verbose ) {
    va_start (args, format);
    fprintf(fp_err, "WARNING : ");
    vfprintf(fp_err, format, args);
    va_end (args);
    fprintf(fp_err, "\n");
  }
}

