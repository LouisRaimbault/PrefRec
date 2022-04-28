#ifndef PREFREC_INIT_H
#define PREFREC_INIT_H


#include <fstream>
#include <string>
#include <iostream>
#include <unordered_map>
#include <vector>

extern const size_t SUL ;
extern const size_t SI  ;
extern const size_t SPN ;
extern const size_t SRN ;


struct pnodes 
{
int sup ;
int litem ;
pnodes * brother;
pnodes * son ;
};

struct rnodes
{   short int indic;
    int sup ;
    uint64_t * tab;
};


void Init_data (char * pathfile, int & nrows, int & maxul ,int & nvar, uint64_t **& Bitdata, std::vector<std::string> & varnames, char deli, double d_sup_init ); 
void sortoutdata (uint64_t** Bdata, int * sumvec, int nb_elem, std::vector<std::string>&  varnames, char ordre);
int init_prefixtree (uint64_t ** Bitdata, uint64_t **& Bdata,std::vector<std::string>&  varnames, int support, int nvar, int maxul, char ordre );
void erase_freq (pnodes* Tree);
void extract_and_erase_freq (pnodes* Tree, std::string str, std::string *  listenom, double nb_freq , int size ,std::ofstream  & flux_set );
void extract_and_erase_freq (pnodes* Tree, std::string str, std::string *  listenom, double nb_freq, int size, std::ofstream  & flux_set, std::ofstream  & flux_rules );


#endif
