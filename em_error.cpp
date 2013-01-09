/*!
 * \author Yu, Kwangmin <yukwangmin@gmail.com> 
 * \date Mon Jul. 13 2012
 * 
 *      
 */


#include <stdio.h>
#include <ctime>
#include <cstring>
#include <iostream>

#include "em_error.h"


using namespace std;

/*
void EM_logger(const ostream &stream, const char filename[], const char functionname[], size_t linenumber, int processID, char contents[]) {
  struct tm *timeinfo;
  time_t t = time(NULL);;
  timeinfo = localtime(&t);
  cout << "[" << asctime(timeinfo) << "] " << filename << " : " << functionname << "() at line " << linenumber << " : [PID:" << processID << "]: " << (contents) << endl;
}
*/

void EM_logger(const ostream &stream, const char filename[], const char functionname[], size_t linenumber, int processID, const string& contents) {
  //sph_logger(filename, functionname, linenumber, processID, contents.c_str());
  //struct tm *timeinfo;
  time_t t = time(NULL);
  #ifndef _MSC_VER
  char *str = ctime(&t);
  #else
  char str[256];
  ctime_s(str, 256, &t);
  #endif
  size_t len = strlen(str);
  str[len - 1] = '\0';

  cout << "[" << str << "] " << filename << " : " << functionname << "() at line " << linenumber << " : [PID:" << processID << "]: " << (contents) << endl;
}