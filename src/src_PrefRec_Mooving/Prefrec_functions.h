#ifndef PREFREC_FUNCTIONS_H
#define PREFREC_FUNCTIONS_H
#include <chrono>
#include "Prefrec_init.h"

extern int minsup;
extern uint64_t * PT_ULTAB;
extern int nb_freq;


void creabitfield ( int *ind,int* mod,uint64_t * ptul, uint64_t * nlgtabl,  int nbca);
void creabitfieldr ( int *ind,int * mod,uint64_t * ptul, uint64_t * nlgtabl,  int nbdases,  int &reste);
void init_Rnodes (pnodes * Tree , int * ngroot , int* tindices ,int * mod, int& nbcases,  int &reste);
void depthwalk (pnodes * Tree, uint64_t * cand, rnodes* tabrnodes);
double Bitprefrec (uint64_t** Bitdata, std::vector<std::string>& varnames, int maxul, int nvar, char ordre, pnodes ** ret_alpha, rnodes ** rnodes_alpha ,int nb_var_moving, int nb_iter, char order,int & nadd, int & ndel);
void rootwalk (pnodes * root );


#endif
