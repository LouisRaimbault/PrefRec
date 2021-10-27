#ifndef PREFIX_TREE_SET_H
#define PREFIX_TREE_SET_H
#include <fstream>
#include <string>
#include <bitset>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <iostream>
#include <bitset>



void Init_data (char * pathfile, int & nrows, int & maxul ,int & nvar, uint64_t **& Bitdata, std::vector<std::string> & varnames, char deli );


struct pnodes 
{
pnodes ( int sup_, int litem_, pnodes* brother_ );
int sup = 0;
unsigned int litem =0;
pnodes * brother = NULL;
pnodes * son = NULL;
};

struct rnodes
{   rnodes (){}
    rnodes ( int sup_,bool indic_, uint64_t* tab_);
     int sup = 0;
    bool indic =0;
    uint64_t * tab = NULL;
};

void sortoutdata (uint64_t** BinaryData, int * sumvec, int nb_elem, std::vector<std::string> &varnames, char ordre);
int init_prefixtree (uint64_t ** Bitdata, uint64_t **& Bdata,std::vector<std::string>& varnames, int support, int nvar, int maxul, char ordre );
void creabitfield ( int *ind,int* mod,uint64_t * ptul, uint64_t * nlgtabl,  int nbc);
void creabitfieldr (  int *ind, int * mod, uint64_t * ptul, uint64_t * nlgtabl,  int nbdases,  int &reste);
void init_Rnodes (rnodes* outrnodes, int* vecsom,  int* tindices ,int * mod,uint64_t** Bitdata, int &nbleft, int& nbcases, int Mins,  int &reste);
void depthwalk (pnodes * Tree, const __int128 * cand, rnodes* tabrnodes , int nbcases);
void rootwalk (pnodes & root ,rnodes * rootnodes,  int nbcases, int nbleft );
pnodes Bitprefrec (uint64_t** Bitdata, std::vector<std::string> &varnames, int maxul, int nvar, char ordre);
void extract_and_erase_freq (pnodes*& Tree);
#endif
