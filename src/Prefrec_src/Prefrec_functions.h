#ifndef PREFREC_FUNCTIONS_H
#define PREFREC_FUNCTIONS_H
#include <chrono>
#include "Prefrec_init.h"

extern int minsup;
extern uint64_t * PT_ULTAB;
extern int nb_freq;


void creabitfield ( int *ind,int* mod,uint64_t * ptul, uint64_t * nlgtabl,  int nbca);
void creabitfieldr ( int *ind,int * mod,uint64_t * ptul, uint64_t * nlgtabl,  int nbdases,  int &reste);
void init_Rnodes (rnodes* outrnodes, int* vecsom,  int* tindices ,int * mod,uint64_t** Bitdata, int &nbleft, int& nbcases, int Mins,  int &reste);
void depthwalk (pnodes * Tree, uint64_t * cand, rnodes* tabrnodes);
void Bitprefrec (uint64_t** Bitdata, std::vector<std::string>& varnames, int maxul, int nvar, char ordre, pnodes *& root);
void rootwalk (pnodes & root ,rnodes * rootnodes,  int nbcases, int nbleft );


#endif
