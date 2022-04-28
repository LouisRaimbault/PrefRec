#include "Prefrules_init.h"


pnodesr * clist = NULL;
pnodesr * ctree = NULL;
pnodesr * cfath = NULL;

pnodesr::pnodesr (double sup_, int litem_, int size_, int is_target_ ) : sup(sup_), litem(litem_), size(size_), is_target(is_target_){};
rules::rules (double conf_, std::string  ante_, std::string cons_, std::string set_) : conf(conf_), ante(ante_), cons(cons_), set(set_) {};


int  double_equals (double a , double b)
{
  return fabs(a - b) < 0.000001;
}

void timestamp()
{
  time_t rawtime;
  struct tm * timeinfo;
  char buffer [30];
  time (&rawtime);
  timeinfo = localtime (&rawtime);

  strftime (buffer,80,"%F. %X : ",timeinfo);
  std::cout << buffer;
}


void init_data (char * pathfile_info, char * pathfile_var, int & nrows, int & nvar, std::vector<std::string> & varnames ,std::vector<double>&  supvalue, std::vector<int>& sizevalue, std::vector<int> & litemvalue )
{ 
  
  long taille = 0; int i = 0; int n = 0; int p = 0; int t = 0;
  char c ; int sit = 0; std::string st="";
  char * tabchar_info = NULL;
  char * tabchar_var = NULL;
  FILE *Fichier = fopen(pathfile_info,"r");
  nrows = 0; 

  if (Fichier != NULL)
    { 
      while (!feof(Fichier))  
       {c = getc(Fichier);
        taille ++;
        if (c == '\n') nrows++; }                      

    }
  fclose(Fichier);
  tabchar_info = (char*)malloc(taille * sizeof(char));
  Fichier = fopen(pathfile_info,"r");
  if ( Fichier != NULL && tabchar_info  != NULL)
    { 
     while (i < taille)
      { c = getc(Fichier);
        tabchar_info[i++]=c;}
      }
  else perror ("\n\n problem in file ");
  fclose(Fichier);

  std::vector<std::string> info_for_tree (nrows);
      
  st.reserve(1000);

  for (t =0; t < nrows; t++, sit++ )
    { while (tabchar_info[sit] != '\n') {st += tabchar_info[sit++];}
          info_for_tree[t]= st;
          st = "";
    }
    free(tabchar_info);

    supvalue = std::vector<double> (nrows);
    sizevalue = std::vector<int> (nrows);
    litemvalue = std::vector<int> (nrows);
    taille = 0; i = 0; sit = 0;
    st.erase(st.begin(),st.end());

    p = 0;
    for (std::string &stk : info_for_tree)
        {
          for (char ch:stk)
            {
                if (ch != '\t') {st+=ch;}
                  else { 
                            if ( n==0 ) {supvalue[p] = std::stod(st); n=1; st.erase(st.begin(),st.end());}
                            else { sizevalue[p] = std::stoi(st); n = 0; st.erase(st.begin(),st.end());}

                        }   
            } 
         litemvalue[p] = std::stoi(st);
         st.erase(st.begin(),st.end());
         p++;    
        }
  

  st.erase(st.begin(),st.end());
  FILE *Fichier_two = fopen(pathfile_var,"r"); 

  if (Fichier_two != NULL)
    { 
      while (!feof(Fichier_two))  
       {c = getc(Fichier_two);
        taille ++;
        if (c == '\n') nvar++; }                      

    }
  fclose(Fichier_two);
  varnames.reserve(nvar);
  tabchar_var = (char*)malloc(taille * sizeof(char));
  Fichier_two = fopen(pathfile_var,"r");
  if ( Fichier_two != NULL && tabchar_var  != NULL)
    { 
     while (i < taille)
      { c = getc(Fichier_two);
        tabchar_var[i++]=c;}
      }
  else perror ("\n\n problem in file ");
  fclose(Fichier_two); 

    for (t =0; t < nvar; t++, sit++ )
    { while (tabchar_var[sit] != '\n') {st += tabchar_var[sit++];}
          varnames.push_back(st);
          st = "";
    }       

free (tabchar_var);

}

void set_coeff_fc (fction_selec * tab_fcs, int * tab_map , std::vector<std::string> & complete_coeff ,int nb_coeff)
{ 
  for (int i = 0; i < nb_coeff;i++)
  { 
    switch (tab_map[i])
    {
    case 0:
      tab_fcs[i].name_fction = complete_coeff[0];
      tab_fcs[i].pfunction = &conftwo;
      break;

    case 1:
      tab_fcs[i].name_fction = complete_coeff[1];
      tab_fcs[i].pfunction = &confone;
      break;

    case 2:
      tab_fcs[i].name_fction = complete_coeff[2];
      tab_fcs[i].pfunction = &powertwo;
      break;
     
      
    case 3:
      tab_fcs[i].name_fction = complete_coeff[3];
      tab_fcs[i].pfunction = &powerone;
      break;

      
    case 4:
      tab_fcs[i].name_fction = complete_coeff[4];
      tab_fcs[i].pfunction = &cov;
      break;

      
    case 5:
      tab_fcs[i].name_fction = complete_coeff[5];
      tab_fcs[i].pfunction = &corr;
      break;

      
    case 6:
      tab_fcs[i].name_fction = complete_coeff[6];
      tab_fcs[i].pfunction = &kappa;
      break;

      
    case 7:
      tab_fcs[i].name_fction = complete_coeff[7];
      tab_fcs[i].pfunction = &maxwellP;
      break;
      
    default :
      tab_fcs[i].name_fction = complete_coeff[0];
      tab_fcs[i].pfunction = &conftwo;      
    }
  }    
}  

int check_target (std::string &set, std::vector<std::string> &targets)
{
  
  std::string stemp = "";
  for (char c: set)
  {
    if (c ==',') 
    {  
      stemp+= c;
      for (std::string & st : targets)
      { 
        if (stemp == st) {return 1;}
      }
      stemp = ""; continue;
    }
    stemp += c;
  }

  return 0;
}


void set_fction_select (fction_selec * tab_fcs, fction_selec * tab_fcs_choice ,std::vector<std::string> & coeff, std::vector<std::string> & coeff_choice ,std::vector<double> & val_coeff)
  { std::vector<std::string> complete_coeff = {"Conf2","Conf1","Power2","Power1","Cov","Corr","Kappa","Maxp"};
    std::unordered_map<std::string,int> map_abuse;
    int i;
    int nb_coeff = coeff.size();
    int nb_coeff_choice = coeff_choice.size();
    for (i = 0; i < complete_coeff.size();i++) {map_abuse[complete_coeff[i]] = i; }

    int tab_sok [nb_coeff];
    int tab_schoice [nb_coeff_choice];
    for ( i = 0; i < nb_coeff;i++)
    { 
      tab_sok [i]= map_abuse[coeff[i]];
      tab_fcs [i].treshold = val_coeff[i];

    }
    
    for ( i = 0; i < nb_coeff_choice;i++)
    { 

      tab_schoice [i]= map_abuse[coeff_choice[i]];
    }
    set_coeff_fc(tab_fcs,tab_sok,complete_coeff,nb_coeff);
    set_coeff_fc(tab_fcs_choice,tab_schoice,complete_coeff,nb_coeff_choice);
  }


void set_target_indices ( bool * ind_target, std::vector<int> & litem_glob , std::vector<int> & size_glob, std::vector<int> & litem_target, int & nb_it)
{ 
  int i,j,k;
  for(i = 0; i < nb_it;i++){ind_target[i]=0;}
  nb_it = 0;
  int v = 1;
  int nb_g = litem_glob.size();
  int nb_t = litem_target.size();
  litem_glob.push_back(-1);
  int cursize = 1;
  int curlitem = 0;
  int tab_h [nb_g];
  for (i = 0; i < nb_g;)
  {curlitem = litem_glob[i];
    v=1;
    for (j =0; j <  nb_t;j++)
    { 
      if (litem_target[j]==curlitem)
      { cursize = size_glob[i];
       
        tab_h[nb_it++] = i;
        v=0;
        i++;
        while (size_glob[i] > cursize)
        {
          tab_h[nb_it++] = i;
          i++;
        }
        
        break;
      }
    }
    if(v!=0) i++;
  }
  for ( i =0; i < nb_it;i++) {ind_target[tab_h[i]] = 1;}
}


void set_target_labels (std::vector<std::string> & nameliste, std::vector<std::string> & targets, std::vector<int> & label_targets)
{ int r,t;
  int name_sl = nameliste.size();
  int nb_target = targets.size();
  
  for (r = 0; r < nb_target;r++)
  { std::string tgt = targets[r];
    
    for ( t = 0;t < name_sl;t++ )
    {
      if (nameliste[t] == tgt )  {label_targets.push_back(t);break;}
    }
  }
}

int give_size_set (std::string & str )
{ int s = 0;
  for (char c : str) {if (c ==',') s++; }
  return s;
}

int check_rules (fction_selec * tab_fction, int nb_fction, std::list<rules> & ruleliste, std::unordered_map<std::string,double> & map_sup)
{  
    
    int i,valid;
    int nb_valid = 0;
    double freq_ante,freq_cons,freq_set;
    std::list<rules>::iterator ite = ruleliste.end();
    std::list<rules>::iterator itb = ruleliste.begin();

     
   for (itb = ruleliste.begin(); itb!=ite;itb++)
   { valid = 1;
     
     freq_ante = map_sup[itb->ante]; freq_cons = map_sup[itb->cons]; freq_set = map_sup[itb->set];
     
      for (i = 0; i < nb_fction;i++)
      { 
          if (tab_fction[i].pfunction(freq_ante,freq_cons,freq_set)<tab_fction[i].treshold){ if (double_equals(tab_fction[i].pfunction(freq_ante,freq_cons,freq_set),tab_fction[i].treshold)) {continue;}  valid = 0; break;}
        }
      if (valid) {itb->valid=1; nb_valid++;}
    }
    return nb_valid;
} 



void fill_rules (fction_selec * tab_fction, int nb_fction, std::list<rules> & ruleliste, std::unordered_map<std::string,double> & map_sup, std::vector<std::vector<double>>&  Mat_coeff,std::vector<std::vector<int>>& Mat_size, std::vector<std::vector<std::string>> & Mat_names, bool ok_sup, bool ok_size )
{ int i=0;int sit = 0; int y = 0;

  double freq_ante,freq_cons,freq_set;
  std::list<rules>::iterator ite = ruleliste.end();
  std::list<rules>::iterator itb = ruleliste.begin();

  for (itb = ruleliste.begin(); itb!=ite;itb++)
  { if(!itb->valid) continue;
    y = 0;
    freq_ante = map_sup[itb->ante]; freq_cons = map_sup[itb->cons]; freq_set = map_sup[itb->set];
    Mat_names[0][sit] = itb->ante; Mat_names[1][sit] = itb->cons; Mat_names[2][sit] = itb->set;
    if(ok_size) {Mat_size[0][sit] = give_size_set(itb->ante); Mat_size[1][sit] = give_size_set(itb->cons); Mat_size[2][sit] = give_size_set(itb->set);}
    if(ok_sup) {Mat_coeff[0][sit] = freq_ante; Mat_coeff[1][sit]=freq_cons; Mat_coeff[2][sit]=freq_set; i = 3;}
    for(i; i < nb_fction;i++)
      { 
        
        Mat_coeff[i][sit] = tab_fction[y].pfunction(freq_ante,freq_cons,freq_set);
        y++;        
      }

    
    sit++;  
   
  }
  
} 


void gen_tree (pnodesr * curlist, pnodesr * curtree, pnodesr * curfather,pnodesr ** listepnodesr, int  & sit, int & stop)
{ if (sit != stop)
{
  if (curlist->size > curtree->size)
  { curlist->father = curtree; 
    curtree->son = curlist;
    sit++;
    gen_tree (listepnodesr[sit],curlist,curtree,listepnodesr,sit,stop);
  }
  if (curlist->size == curtree->size)
  {   curlist->father = curfather;
    curtree->brother = curlist;
    sit++;
    gen_tree (listepnodesr[sit],curlist,curfather,listepnodesr,sit,stop);
    
  }
  
  if (curlist->size < curtree->size)
  {
    gen_tree (curlist, curtree->father,curtree->father->father,listepnodesr,sit,stop);
  }
}
else {clist = curlist; ctree = curtree; cfath = curfather;}   

}


void rebuild_tree (std::vector<double> & supvalue, std::vector<int> & sizevalue, std::vector<int> & litemvalue, bool * target_tab ,pnodesr *& rootrules, pnodesr** tabpnodes ,int nb_freq )
{
  int sit = 0; int nbg = 0; int j; 
  int k = nb_freq/10000;
  int reste_k = nb_freq-10000*k;
  int ct = 1;


  for (j = 0; j < nb_freq;j++)
    { 
      tabpnodes[j] = new pnodesr(supvalue[j],litemvalue[j],sizevalue[j],target_tab[j]);
    }
  
  clist = tabpnodes[0];
  ctree = rootrules;
  cfath = rootrules;
  for(j = 0; j< k; j++,ct++)
    {
      nbg = 10000 * ct;
      gen_tree(clist,ctree,cfath,tabpnodes,sit,nbg);
    }
  gen_tree(clist,ctree,cfath,tabpnodes,sit,nb_freq);
}


void erase_tree (pnodesr *& Tree)
{ 
  if (Tree->son) {erase_tree(Tree->son);}
  if (Tree->brother) {erase_tree(Tree->brother);}
  delete Tree;
}