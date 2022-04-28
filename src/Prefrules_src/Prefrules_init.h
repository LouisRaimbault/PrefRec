#ifndef PREFRULES_INIT_H
#define PREFRULES_INIT_H


#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <unordered_map>
#include <iostream>
#include <time.h>
#include "Prefrules_coefficient.h"


struct pnodesr
{
 pnodesr (double sup_, int litem_, int size_, int is_target_ );
  pnodesr * brother = NULL;
  pnodesr * son = NULL;
  pnodesr * father = NULL; 
  double sup = 0;
  bool is_target = 0;
  int size =0;
  int litem =0;
};


struct rules
  { 
    rules (double conf_, std::string  ante_, std::string cons_, std::string set_);
    double conf =0;
    std::string ante ;
    std::string cons ;
    std::string set ;
    bool valid = 0;
  };

struct fction_selec
{
  double (*pfunction)(double freq_ante, double freq_complem, double freq_set );
  double treshold = 0;
  double max = 0;
  std::string name_fction = "";  
};


extern pnodesr * clist;
extern pnodesr * ctree;
extern pnodesr * cfath;
 

void timestamp();
void init_data (char * pathfile_info, char * pathfile_var, int & nrows, int & nvar, std::vector<std::string> & varnames ,std::vector<double>&  supvalue, std::vector<int> & sizevalue, std::vector<int> & litemvalue );
int st_toi (std::string coeff);
int double_equals(double a , double b);
void set_coeff_fc (fction_selec * tab_fcs, int * tab_map , std::vector<std::string> & complete_coeff ,int nb_coeff);
void set_target_indices ( bool * ind_target, std::vector<int> & litem_glob , std::vector<int> & size_glob, std::vector<int> & litem_target, int & nb_it);
int check_target (std::string &set, std::vector<std::string> &targets);
void set_target_labels (std::vector<std::string> & nameliste, std::vector<std::string> & targets, std::vector<int> & label_targets);
void set_fction_select (fction_selec * tab_fcs, fction_selec * tab_fcs_choice ,std::vector<std::string> & coeff, std::vector<std::string> & coeff_choice ,std::vector<double> & val_coeff);
void erase_tree (pnodesr *& Tree);
void rebuild_tree (std::vector<double> & supvalue, std::vector<int> & sizevalue, std::vector<int> & litemvalue, bool * target_tab ,pnodesr *& rootrules, pnodesr** tabpnodes ,int nb_freq );
void gen_tree (pnodesr * curlist, pnodesr * curtree, pnodesr * curfather,pnodesr ** listepnodesr, int  & sit, int & stop);
int check_rules (fction_selec * tab_fction, int nb_fction, std::list<rules> & ruleliste, std::unordered_map<std::string,double> & map_sup);
void fill_rules (fction_selec * tab_fction, int nb_fction, std::list<rules> & ruleliste, std::unordered_map<std::string,double> & map_sup, std::vector<std::vector<double>>&  Mat_coeff,std::vector<std::vector<int>>& Mat_size, std::vector<std::vector<std::string>> & Mat_names, bool ok_sup, bool ok_size );





#endif