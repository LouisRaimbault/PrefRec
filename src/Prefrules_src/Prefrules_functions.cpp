

#include "Prefrules_functions.h"


void allrules  (std::string & stree, std::string& cntree , std::string * tabconseqalpha, int size , int good, std::unordered_map<std::string,double>& map_double, double cursup, int cursize, double Minconf,std::list<rules>& ruliste, std::string & strset )
{ 
  
  std::string trstree = stree;
  std::string cnstree = cntree;
  int ngr = 0;
  int tail = size - good;
  
  if (tail == 1)
  {  
    trstree.erase(trstree.find(tabconseqalpha[good]),tabconseqalpha[good].size());
    double indic = conftwo(map_double[trstree],map_double[cnstree],cursup);
    if ((indic < Minconf == 0))
    { 
      cnstree +=  tabconseqalpha[good];
      ruliste.emplace_front (rules(indic,trstree,cnstree,strset));
    }
    
    
  }
  else
  { 
    std::string * antecedant = new std::string [tail];
    std::string * consequent  = new std::string [tail];
    std::string * consalpha   = new std::string [tail];
    
    for (int i = good; i < size ; i++)
    {
      
      trstree.erase(trstree.find(tabconseqalpha[i]),tabconseqalpha[i].size());
      double indic = conftwo(map_double[trstree],map_double[cnstree],cursup);
      if ( (indic < Minconf) == 0)
      { 
        cnstree +=  tabconseqalpha[i];
        antecedant[ngr] = trstree;
        consequent[ngr] = cnstree;
        consalpha[ngr++]= tabconseqalpha[i];
        ruliste.emplace_front (rules (indic,trstree,cnstree,strset));
      }
      trstree = stree; 
      cnstree = cntree;
    }
    
    if (ngr >1 && cursize> 2) 
    {  
      int ngt = ngr -1;
      for (int j = 0; j < ngt;j++)
      { 
        allrules(antecedant[j],consequent[j],consalpha,ngr,j+1,map_double,cursup,cursize-1,Minconf, ruliste,strset);
      }
    }
    
    delete [] antecedant;
    delete [] consequent;
    delete [] consalpha;
    
  }
  
  
}


void genrules (pnodesr & tree ,int size ,std::string ** tabname,std::string & str, std::string * varnames, std::unordered_map<std::string,double>& map_double, double& Minconf, std::list<rules>& ruliste)
{ 
  std::string strtemp = str + varnames[tree.litem];
  double cursup = tree.sup;
  map_double[strtemp] = tree.sup;
  tabname[size-1] = &varnames[tree.litem];
  if(tree.is_target){
  
  std::string * antecedant = new std::string [size];
  std::string * consequent  = new std::string [size];
  std::string * consalpha   =  new std::string [size];
  int nrg = 0;     
  std::string tmpstr ;
  tmpstr.reserve(100);
  tmpstr +=  *(tabname[1]);
  for (int i = 2; i < size;i++)
  { tmpstr += *(tabname[i]);}
  double indic = conftwo(map_double[tmpstr],map_double[*tabname[0]],cursup);
  if ( (indic < Minconf) == 0 ) {
    ruliste.emplace_front(rules(indic,tmpstr,*tabname[0],strtemp));
    antecedant[nrg] = tmpstr;
    consalpha[nrg] = *tabname[0];
    consequent [nrg++] = *tabname[0];}
  
  tmpstr.erase(tmpstr.begin(),tmpstr.end());
  tmpstr += *(tabname[0]);
  
  size_t st = tabname[0]->size();
  int l = 0;
  for (int j = 1; j < size; j++)
  {
    for (l =1; l < size ; l++)
    {
      if (l != j) tmpstr += *(tabname[l]);     
    }
    double indic = conftwo(map_double[tmpstr], map_double[*tabname[j]],cursup);
    if ((indic < Minconf) == 0)  
    {
      ruliste.emplace_front(rules(indic,tmpstr,*tabname[j],strtemp));
      antecedant[nrg] = tmpstr;
      consalpha [nrg] = *tabname[j];
      consequent[nrg++] = *tabname[j];
    } 
    tmpstr.erase(tmpstr.begin()+st,tmpstr.end());       
  }
  
  if (nrg >1 && size > 2) 
  { int ngt = nrg -1;
    for (int m = 0; m < ngt;m++)
    { 
      allrules(antecedant[m],consequent[m],consalpha,nrg,m+1,map_double,cursup,size-1,Minconf, ruliste,strtemp);
    }
  }
  delete [] antecedant;
  delete [] consequent;
  delete [] consalpha; 
  }
  
  if (tree.son) 
    {
     genrules(*tree.son,size+1,tabname,strtemp,varnames,map_double,Minconf,ruliste);
     }
  
  if (tree.brother)
    { 
      genrules(*tree.brother,size,tabname,str,varnames,map_double,Minconf,ruliste);
    }

} 


void Genrules_root (pnodesr & roots, int size, std::string ** tabname, std::string &str, std::string * varnames, std::unordered_map<std::string,double>& map_double,double & Minconf, std::list<rules>& ruliste)
{ 
  std::string strtemp = str + varnames[roots.litem];
  double cursup = roots.sup;
  map_double[strtemp] = cursup;
  tabname[0] = &varnames[roots.litem];
  if (roots.son)
  { genrules (*roots.son,size+1,tabname,strtemp,varnames,map_double, Minconf,ruliste);}
  if (roots.brother) {  Genrules_root(*roots.brother,size,tabname,str,varnames,map_double, Minconf,ruliste);}    
}

