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


struct dimentions    // les dimentions de la matrice
{ public:
  dimentions(int nb_elem_, int nb_lignes_);
  int nb_elem = 0;
  int nb_lignes = 0;
  int nb_var = 0; 
};

int sumbool (bool* tab, int n);

dimentions longueurfichier (char * cheminfichier); //détermine la longueur d'un fichier 

void remplirfichier (char *cheminfichier, char* chaine, dimentions& dimtableau); // prends une chaine de character ayant la bonne chaine 
void getmap  (std::vector<std::string> transac, std::unordered_map<std::string,int>& maptr,std::vector<std::string> &listename ,char deli);
void transformintotransac (char* chainecara, std::string* outtransac,dimentions& dimtab);
void transactiontoBitmax( bool ** mytab,std::unordered_map<std::string,int> &mappy,std::vector<std::string> transac, std::vector<std::string> &listename ,char deli,int nb_individus);
void tri_tableau (bool** dataframept, std::vector<int> mytab, int nb_elem, std::vector<std::string> &listenoms);

int init_prefixtree (bool** mydataframe,std::vector<std::string>& listenoms, int support, int nbvar, int nbin );
void Bitmaxtranspos (bool ** Bitranspose,bool ** Bitmax, int nbc, int nbl);

struct pnodes 
{
pnodes ( int sup_, unsigned int litem_, pnodes* brother_ );
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

void rootsumtest ( int * outsum, int* indyces, bool ** TRBitmax , int supcand, int nbleftrn);
void getlongfull (uint64_t * nod, uint64_t * cand, uint64_t* nouvroot, int & nbelem);
void transformintolong ( int *ind,bool * ptshort, uint64_t * nlgtabl,  int nbc);
void transformintolongreste ( int *ind,bool * ptshort, uint64_t * nlgtabl,  int nbdases,  int &reste);
void createvecrnodes (rnodes* outrnodes, int* vecsom,  int* tindices ,bool** Bitmax, int &nbleft, int& nbcases, int Mins,  int &reste);
void depthwalk (pnodes & Tree, uint64_t * cand, rnodes* tabrnodes ,  int nbcases);
void rootwalk (pnodes & root ,rnodes * rootnodes, unsigned int nbcases, int nbleft );
pnodes Bitprefrec (bool** Bitmax, std::vector<std::string> &varnames, double relativeSup, int nbind, int nbvar);
void extract_and_erase_freq (pnodes*& Tree, std::string str, std::string *  listenom, std::string * namevalue, int * supvalue, int* sizevalue, int* litemvalue , int size ,int& situeur );


#endif