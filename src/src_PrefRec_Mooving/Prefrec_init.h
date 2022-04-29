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


struct rnodes;

struct pnodes 
{
int sup ;
rnodes * link_rn;
pnodes * brother;
pnodes * son ;
};

struct rnodes
{
    int sup_data ;
    int sup;
    uint64_t * tab;
    rnodes * brother;
    pnodes * cor_bsb; 
    uint64_t * tab_data;
    std::string * var;
    int lab;
};


void Init_data (char * pathfile, int & nrows, int & maxul ,int & nvar, uint64_t **& Bitdata, std::vector<std::string> & varnames, char deli ); 
void sortoutdata (uint64_t** Bdata, int * sumvec, int nb_elem, std::vector<std::string>&  varnames, char order);
int init_prefixtree (uint64_t ** Bitdata, uint64_t **& Bdata,std::vector<std::string>&  varnames, int support, int nvar, int maxul, char order);
void erase_freq (pnodes* Tree);
void extract_and_erase_freq (pnodes* Tree, std::string str, std::string *  listenom, double nb_freq , int size ,std::ofstream  & flux_set );
void extract_and_erase_freq (pnodes* Tree, std::string str, std::string *  listenom, double nb_freq, int size, std::ofstream  & flux_set, std::ofstream  & flux_rules );
void set_rnodes_in_order (rnodes * root_rnodes, rnodes * n_rnodes);
void erase_condition_father (pnodes * pn, rnodes  * n_rnodes, int * nbdel);
void erase_mobile (pnodes * pn, int * nbdel);
void erase_var (rnodes * rn, int * nbdel);
void erase_nouvar (pnodes * Tree, rnodes * rn ,int * nbdel);


#endif
