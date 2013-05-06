/*
 * Backtrace.h
 *
 *  Created on: Apr 24, 2013
 *      Author: martin
 */

#ifndef BACKTRACE_H_
#define BACKTRACE_H_

#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <exception>
#include <cxxabi.h>
#include <string.h>

class Backtrace
{
private:
  void *array[10];
  size_t size;

  std::string
  demangle(const char* symbol) const
  {

    int status;
    char* demangled;

    //allocate mem

    char *tmp = (char *) malloc(strlen(symbol) * sizeof(char));

    if (tmp == NULL)
      {
        return symbol;
      }

    //first, try to demangle a c++ name

    if (1 == sscanf(symbol, "%*[^(]%*[^_]%[^)+]", tmp))
      {
        demangled = abi::__cxa_demangle(tmp, NULL, NULL, &status);
        if (status == 0)
          {
            std::string result(demangled);
            free(demangled);
            return result;
          }
      }

    return symbol;
  }

public:

  Backtrace()
  {
    size = backtrace(array, 10);
  }

  void
  print_trace() const
  {
    // print out all the frames to stderr

                char** symbols = backtrace_symbols(array, size);

                for (unsigned int i = 1; i < size; ++i) {
                        printf("%i) %s\n\n", i, demangle(symbols[i]).c_str());
                }

                free(symbols);
  }
};

#endif /* BACKTRACE_H_ */
