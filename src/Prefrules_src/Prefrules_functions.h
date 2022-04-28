#ifndef PREFRULES_FUNCTIONS_H
#define PREFRULES_FUNCTIONS_H

#include "Prefrules_init.h"


void allrules  (std::string & stree, std::string& cntree , std::string * tabconseqalpha, int size , int good, std::unordered_map<std::string,double>& mappy, double cursup, int cursize, double Minconf,std::list<rules>& ruliste, std::string & strset);

void genrules (pnodesr & tree ,int size ,std::string ** tabname,std::string & str, std::string * varnames, std::unordered_map<std::string,double>& mappy, double& Minconf, std::list<rules>& ruliste );

void Genrules_root (pnodesr & roots, int size, std::string ** tabname, std::string &str, std::string * varnames, std::unordered_map<std::string,double>& map_double,double & Minconf, std::list<rules>& ruliste);




#endif
