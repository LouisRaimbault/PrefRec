#include <fstream>
#include <string>
#include <bitset>
#include <vector>
#include <iostream>
#include <list>
#include <unordered_map>
#include <iostream>
#include <bitset>
#include <iostream>
#include <math.h>

struct dimentions    // les dimentions de la matrice
{ public:
  dimentions(int nb_elem_, int nb_lignes_);
  int nb_elem = 0;
  int nb_lignes = 0;
  int nb_var = 0; 
};

struct pnodesr
{
 pnodesr (double sup_, int litem_, int size_ );
 double sup = 0;
 int litem =0;
 pnodesr * brother = NULL;
 pnodesr * son = NULL;
 pnodesr * father = NULL; 
 int size =0;
};


struct rules 
  { rules (double conf_, std::string  ante_, std::string comp_, std::string set_);
    double conf =0;
    std::string ante ;
    std::string comp ;
    std::string set ;
  };


struct fctionpt 
  {
  double (*mafonc)(double & freq_ante, double & freq_complem, double & freq_set );
  };

dimentions longueurfichier (char * cheminfichier); //détermine la longueur d'un fichier 

void remplirfichier (char *cheminfichier, char* chaine, dimentions& dimtableau); // prends une chaine de character ayant la bonne chaine 

void transformintotransac (char* chainecara, std::string * pstring ,dimentions& dimtab);

void getabfull (std::vector<std::string>& tabstr, double * RelativeSupvalue, int *  sizevalue, int *  litemvalue);

double conftwo (double & freq_ante, double & freq_complem, double & freq_set );

double powertwo (double & freq_ante, double & freq_complem, double & freq_set );

double confone (double & freq_ante, double & freq_complem, double & freq_set );

double powerone (double & freq_ante, double & freq_complem, double & freq_set );

double cov (double & freq_ante, double & freq_complem, double & freq_set );

double corr (double & freq_ante, double & freq_complem, double & freq_set);

double kappa (double & freq_ante, double & freq_complem, double & freq_set);

double MaxwellP (double & freq_ante, double & freq_complem, double & freq_set);

void gen_tree (pnodesr * curlist, pnodesr * curtree, pnodesr * curfather,pnodesr ** listepnodesr, int  & sit, int & stop);

void gen_listetree (std::vector<double> & supvalue, std::vector<int> & litemvalue, std::vector<int> & sizevalue, int nbfreq, pnodesr** tabpnodes);

void allrules  (std::string & stree, std::string& cntree , std::string * tabconseqalpha, int size , int good, std::unordered_map<std::string,double>& mappy, double cursup, int cursize, double Minconf,std::list<rules>& ruliste, std::string & strset);

void genrules (pnodesr & tree ,int size ,std::string ** tabname,std::string & str, std::string * varnames, std::unordered_map<std::string,double>& mappy, double& Minconf, std::list<rules>& ruliste );

void Genrulesroot (pnodesr & roots, int size, std::string ** tabname, std::string &str, std::string * varnames, std::unordered_map<std::string,double>& mappy,double & Minconf, std::list<rules>& ruliste);



