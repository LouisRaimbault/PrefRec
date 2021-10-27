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

void init_data (char * pathfile_info, char * pathfile_var, int & nrows, int & nvar, std::vector<std::string> & varnames ,double *& supvalue, int *& sizevalue, int *&litemvalue );

double conftwo (double & freq_ante, double & freq_complem, double & freq_set );

double cov (double & freq_ante, double & freq_complem, double & freq_set );

double corr (double & freq_ante, double & freq_complem, double & freq_set);

double kappa (double & freq_ante, double & freq_complem, double & freq_set);

double MaxwellP (double & freq_ante, double & freq_complem, double & freq_set);

void gen_tree (pnodesr * curlist, pnodesr * curtree, pnodesr * curfather,pnodesr ** listepnodesr, int  & sit, int & stop);

void gen_listetree (double *  supvalue, int * litemvalue, int * sizevalue, int nbfreq, pnodesr** tabpnodes);

void allrules  (std::string & stree, std::string& cntree , std::string * tabconseqalpha, int size , int good, std::unordered_map<std::string,double>& mappy, double cursup, int cursize, double Minconf,std::list<rules>& ruliste, std::string & strset);

void genrules (pnodesr & tree ,int size ,std::string ** tabname,std::string & str, std::string * varnames, std::unordered_map<std::string,double>& mappy, double& Minconf, std::list<rules>& ruliste );

void Genrulesroot (pnodesr & roots, int size, std::string ** tabname, std::string &str, std::string * varnames, std::unordered_map<std::string,double>& mappy,double & Minconf, std::list<rules>& ruliste);



