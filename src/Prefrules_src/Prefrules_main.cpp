#include "Prefrules_init.h"
#include "Prefrules_functions.h"
#include "Prefrules_coefficient.h"


int main (int argc , char ** argv)
{

   timestamp(); 
   std::cout << "Checking parameters ... \n"; 
   std::vector<std::string> complete_coeff = {"Conf2","Conf1","Power2","Power1","Cov","Corr","Kappa","Maxp"};
   int n_coeff_pos = 8;
   int ar,i,j,k,n,r,npa,nb_freq;
   int ok_target = 1;
   std::vector<double> supvalue;
   std::vector<int> sizevalue;
   std::vector<int> litemvalue;
   std::vector<std::string> varnames;
   int nvar =0; int nrows = 0;

   char * pathr = (char*)argv[1];
   char * pathv = (char*)argv[2];
   std::vector<std::string> targets;
   int nb_coeff_test = std::atoi((char*)argv[3]);
   std::vector<std::string> coeff_test (nb_coeff_test);
   std::vector<double> value_coeff_test (nb_coeff_test);
   ar = 4;
   for (i=0; i < nb_coeff_test; i++,ar++)
    {
      std::string tr_ct = std::string((char*)argv[ar]);
      std::string tr_c = "";
      for (j = 0; j < tr_ct.size();j++)
        {
          if (tr_ct[j] == '=') {coeff_test[i] = tr_c; j++; tr_c = ""; break;}
          tr_c += tr_ct[j];
        } 

      for (j; j < tr_ct.size();j++)
        {
          tr_c += tr_ct[j];
        }

      value_coeff_test[i] = std::stod(tr_c); 
    }
    if (coeff_test[0]!= "Conf2") {std::cout << "first coefficient condition must be Conf2 \n"; return 0;}
    for (i = 0; i < nb_coeff_test;i++) 
      { 
        k = 0;
        for (j = 0; j < n_coeff_pos;j++)
          { 
            if (coeff_test[i] == complete_coeff[j] ) {k=1; break;}
          }
        if (k == 0) {std::cout << "Error : " << coeff_test[i] << " is not a possible coefficient \n"; return 0; }
      }


    int nb_coeff_extract = std::atoi((char*)argv[ar]);
    std::vector<std::string> coeff_extract (nb_coeff_extract);
    ar++;
    for (i = 0;i < nb_coeff_extract;i++,ar++)
      {
        coeff_extract[i] = std::string((char*)argv[ar]);
      }
    if (coeff_extract[0] == "all_indicators") {nb_coeff_extract = 8; coeff_extract = std::vector<std::string> {"Conf2","Conf1","Power2","Power1","Cov","Corr","Kappa","Maxp"};}
    
    std::string target = std::string((char*)argv[ar++]);
    
    if (target == "all_tgts") {ok_target = 0;}
    targets.push_back(target);
    int ok_choice [3];

    for (i = 0; i < 3;i++,ar++)
      {
        ok_choice[i] = std::atoi((char*)argv[ar]);
      }
    timestamp();  
    std::cout << "Frequent itemset data retrieval ... \n";
    init_data(pathr, pathv, nb_freq,nvar, varnames,supvalue,sizevalue,litemvalue);
    std::vector<int> litem_targets;

    bool * ind_target = new bool [nb_freq];
    for (r = 0; r < nb_freq;r++) {ind_target[r] = 1;}

    int nb_target_frequent = nb_freq;
    pnodesr ** tabpnodes = new pnodesr* [nb_freq];
    pnodesr * rootrules =  new pnodesr(0,0,0,0);
    if (ok_target) 
      {  
        set_target_labels(varnames,targets,litem_targets);
        set_target_indices(ind_target,litemvalue,sizevalue,litem_targets,nb_target_frequent);
  
        if (nb_freq == 0)
          {
            std::cout << "No one rule with this target(s) \n";
            return 0;  
          }

      }
  timestamp();
  std::cout << "Constructing Frequent Tree with target :  " << target << " \n";
  rebuild_tree(supvalue,sizevalue,litemvalue,ind_target,rootrules,tabpnodes,nb_freq);
  
  delete [] ind_target;
  
  nb_freq = nb_target_frequent;

  std::unordered_map<std::string,double> map_sup;
  map_sup.reserve(nb_freq);
  fction_selec tab_choice_coeff [nb_coeff_test+1];
  tab_choice_coeff[0].pfunction = &conftwo;
  fction_selec tab_choice_extract [nb_coeff_extract+1];
  
  set_fction_select (tab_choice_coeff,tab_choice_extract,coeff_test,coeff_extract,value_coeff_test);

  
  
  std::vector<std::string> names_size {"Size_Ante","Size_Cons","Size_Set"};
  std::vector<std::string> names_sup {"Sup_Ante","Sup_Cons","Sup_Set"};

  std::list<rules> ruleliste;
  int taille_nameliste = varnames.size();
  std::string * tabname  [taille_nameliste];
  std::string st;
  st.erase(st.begin(),st.end());
  st.reserve (100);
  int compte = 0;
  double minConf = value_coeff_test[0];
  timestamp();
  std::cout << "Search rules with conf2 =  " << minConf << "\n";
  Genrules_root (*rootrules->son,1,tabname,st,&varnames[0],map_sup, minConf ,ruleliste);
  erase_tree(rootrules);
  timestamp();
  std::cout << "Total of Rules :  " << ruleliste.size() << " checking other conditions : ";
  for (i = 1; i < nb_coeff_test;i++)
    {std::cout << coeff_test[i] << " = " << value_coeff_test[i] << " ";}
  std::cout << "\n";
  int n_rules = check_rules(tab_choice_coeff,nb_coeff_test,ruleliste,map_sup);
  timestamp();

  std::cout << "Number of valid Rules " << n_rules << " Rules \n"; 

   if (argc <= ar ) {return 0;}

  timestamp();
  std::cout << "Preparing data to export \n";


  int ok_size = ok_choice[0];
  int ok_sup = ok_choice[1];
  int ok_global_indic = ok_choice[2];
  int total_coeff_double = 0;


  std::vector<std::vector<std::string>> Mat_names (3, std::vector<std::string>(n_rules)); 
  std::vector<std::vector<int>> Mat_size;
  std::vector<std::vector<double>> Mat_coeff;
  
  if (ok_size) {Mat_size = std::vector<std::vector<int>> (3, std::vector<int> (n_rules));}
  if (ok_sup && ok_global_indic )  {total_coeff_double = nb_coeff_extract+5; Mat_coeff = std::vector<std::vector<double>> (nb_coeff_extract+5, std::vector<double> (n_rules));}
  if (ok_sup && !ok_global_indic )  {total_coeff_double = nb_coeff_extract+3; Mat_coeff = std::vector<std::vector<double>> (nb_coeff_extract+3, std::vector<double> (n_rules));}
  if (!ok_sup && ok_global_indic )  {total_coeff_double = nb_coeff_extract+2; Mat_coeff = std::vector<std::vector<double>> (nb_coeff_extract+2, std::vector<double> (n_rules));}
  if (!ok_sup && !ok_global_indic )  {total_coeff_double = nb_coeff_extract;Mat_coeff = std::vector<std::vector<double>> (nb_coeff_extract, std::vector<double> (n_rules));}


  fill_rules(tab_choice_extract,total_coeff_double-3*ok_sup-2*ok_global_indic,ruleliste, map_sup,Mat_coeff,Mat_size,Mat_names,ok_size,ok_sup);

  int nb_frequent_target = 0;
  for (std::unordered_map<std::string,double>::iterator it = map_sup.begin();it != map_sup.end();it++)
  { std::string sttest = it->first;
    nb_frequent_target += check_target(sttest,targets);
  }

  if (ok_global_indic)
    { std::list<rules>::iterator it = ruleliste.begin();
      double * pt_meanc = &Mat_coeff[Mat_coeff.size()-2][0];
      double * pt_strong = &Mat_coeff[Mat_coeff.size()-1][0];
      for (r = 0; r < n_rules; it++)
      { if (!it->valid) continue;
          pt_meanc[r] = (conftwo(map_sup[it->ante],map_sup[it->cons],map_sup[it->set])+confone(map_sup[it->ante],map_sup[it->cons],map_sup[it->set])+powertwo(map_sup[it->ante],map_sup[it->cons],map_sup[it->set])+powerone(map_sup[it->ante],map_sup[it->cons],map_sup[it->set]))/(double)4;
          pt_strong[r] = (double)(1-map_sup[it->ante]-map_sup[it->cons]+2*map_sup[it->set]);
          r++;
        }
    }
  
  
  std::vector<int> target_antecedant;
  std::vector<int> target_consequent;
  
  if (ok_target)
  { 
    
    target_antecedant= std::vector<int> (n_rules);
    target_consequent= std::vector<int> (n_rules);
    std::string * pt_tante = &Mat_names[0][0];
    std::string * pt_tcons = &Mat_names[1][0];
    
    for (r =0; r < n_rules;r++)
      { 
        target_antecedant[r]=check_target(pt_tante[r],targets);
        target_consequent[r]=check_target(pt_tcons[r],targets);
      }
    
  }

 

  std::string pathout = std::string((char*)argv[ar]);
  timestamp();
  std::cout << "writing in file : " << pathout << "\n"; 
  std::ofstream outfluxrules;
  outfluxrules.open(std::string(pathout));
  outfluxrules << "Antecedant\tConsequent\tset";
  if (ok_size) {outfluxrules << "\tSize_Ante\tSize_Cons\tSize_Set";}
  if (ok_sup)  {outfluxrules << "\tSup_Ante\tSup_Cons\tSup_Set";}
  for (i = 0; i < nb_coeff_extract;i++) {outfluxrules << "\t" << coeff_extract[i];}
  if (ok_global_indic) {outfluxrules << "\tMean" << "\tStrong";}
  if (ok_target) {outfluxrules << "\tin_antecedant\tin_consequent";}
  outfluxrules << "\n";


  std::string * pt_antecedant = &Mat_names[0][0];
  std::string * pt_consequent = &Mat_names[1][0];
  std::string * pt_set = &Mat_names[2][0];

  for (i = 0; i < n_rules;i++)
    {
      outfluxrules << pt_antecedant[i] << "\t"<<pt_consequent[i]<<"\t"<<pt_set[i];
      if (ok_size) 
        { 
          for (j = 0; j < 3; j++) {outfluxrules << "\t" << Mat_size[j][i];}   
        }
      for (k = 0; k < total_coeff_double;k++) {outfluxrules << "\t" << Mat_coeff[k][i];}

      if (ok_target) {outfluxrules << "\t"<<target_antecedant[i]<<"\t"<<target_consequent[i];}
      outfluxrules << "\n";
    }
   
  outfluxrules.close();





   return 0;

}