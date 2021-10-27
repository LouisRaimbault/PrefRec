#include "PrefRules.h"



rules::rules (double conf_, std::string  ante_, std::string comp_, std::string set_ ) : conf(conf_), ante(ante_), comp(comp_), set(set_) {};
pnodesr::pnodesr (double sup_, int litem_, int size_ ) : sup(sup_), litem(litem_), size(size_){};



void init_data (char * pathfile_info, char * pathfile_var, int & nrows, int & nvar, std::vector<std::string> & varnames ,double *& supvalue, int *& sizevalue, int *&litemvalue )
{ 
  std::cout << "entree init \n";
  long taille = 0; int i = 0; int n = 0; int p = 0; int t = 0;
  char c ; int sit = 0; std::string st="";
  char * tabchar_info = NULL;
  char * tabchar_var = NULL;
  FILE *Fichier = fopen(pathfile_info,"r"); 

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
    supvalue = (double*)malloc(nrows*sizeof(double));
    sizevalue = (int*)malloc(nrows*sizeof(int));
    litemvalue = (int*)malloc(nrows*sizeof(int));
    taille = 0; i = 0; sit = 0;
    st.erase(st.begin(),st.end());
    for (std::string &stk : info_for_tree)
        {
          for (char ch:stk)
            {
                if (ch != '\t') {st+=ch;}
                  else { 
                            if ( n==0 ) {supvalue[p] = std::stod(st); n=1; st.erase(st.begin(),st.end());}
                            else {sizevalue[p] = std::stoi(st); n = 0; st.erase(st.begin(),st.end());}

                        }   
            }
         litemvalue[p] = std::stoi(st);
         st.erase(st.begin(),st.end());
         p++;    
        }
  std::cout << "partie 2  \n";
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





double conftwo (double & freq_ante, double & freq_complem, double & freq_set )
{
  return freq_set/freq_ante;
  
}


double cov (double & freq_ante, double & freq_complem, double & freq_set )
{
  return freq_set - freq_ante*freq_complem;
}

double corr (double & freq_ante, double & freq_complem, double & freq_set)
{ 
  double cov = freq_set - freq_ante*freq_complem;
  return cov/std::pow(freq_ante*(1-freq_ante)*freq_complem*(1-freq_complem),0.5);
}


double kappa (double & freq_ante, double & freq_complem, double & freq_set)
{ 
  double cov = freq_set - freq_ante*freq_complem;
  return cov/(freq_ante*(1-freq_complem)+freq_complem*(1-freq_ante));
}

double MaxwellP (double & freq_ante, double & freq_complem, double & freq_set)
{
  double cov = freq_set - freq_ante*freq_complem;
  return (2*cov)/(freq_ante*(1-freq_ante)+freq_complem*(1-freq_complem));
}

pnodesr * clist = NULL;
pnodesr * ctree = NULL;
pnodesr * cfath = NULL;

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



void gen_listetree (double *  supvalue, int *  litemvalue, int *  sizevalue, int nbfreq, pnodesr** tabpnodes)
  {
    for (int i =0; i < nbfreq; i++)
      {
        tabpnodes[i] = new pnodesr (supvalue[i],litemvalue[i],sizevalue[i]);
      }
  }




void allrules  (std::string & stree, std::string& cntree , std::string * tabconseqalpha, int size , int good, std::unordered_map<std::string,double>& mappy, double cursup, int cursize, double Minconf,std::list<rules>& ruliste, std::string & strset )
{ 
  
  std::string trstree = stree;
  std::string cnstree = cntree;
  int ngr = 0;
  int tail = size - good;
  
  if (tail == 1)
  {  
    trstree.erase(trstree.find(tabconseqalpha[good]),tabconseqalpha[good].size());
    double indic = conftwo(mappy[trstree],mappy[cnstree],cursup);
    if (indic > Minconf)
    { 
      cnstree +=  tabconseqalpha[good];
      ruliste.emplace_front (rules (indic,trstree,cnstree,strset));
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
      double indic = conftwo(mappy[trstree],mappy[cnstree],cursup);
      if (indic > Minconf)
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
        allrules(antecedant[j],consequent[j],consalpha,ngr,j+1,mappy,cursup,cursize-1,Minconf, ruliste,strset);
      }
    }
    
    delete [] antecedant;
    delete [] consequent;
    delete [] consalpha;
    
  }
  
  
}


void genrules (pnodesr & tree ,int size ,std::string ** tabname,std::string & str, std::string * varnames, std::unordered_map<std::string,double>& mappy, double& Minconf, std::list<rules>& ruliste)
{ 
  std::string strtemp = str + varnames[tree.litem];
  double cursup = tree.sup;
  mappy[strtemp] = tree.sup;
  tabname[size-1] = &varnames[tree.litem];
  std::string * antecedant = new std::string [size];
  std::string * consequent  = new std::string [size];
  std::string * consalpha   =  new std::string [size];
  int nrg = 0;     
  std::string tmpstr ;
  tmpstr.reserve(50);
  tmpstr +=  *(tabname[1]);
  for (int i = 2; i < size;i++)
  { tmpstr += *(tabname[i]);}
  double indic = conftwo(mappy[tmpstr],mappy[*tabname[0]],cursup);
  if (indic > Minconf) {
    ruliste.emplace_front(rules(indic,tmpstr,*tabname[0],strtemp));
    antecedant[nrg] = tmpstr;
    consalpha[nrg] = *tabname[0];
    consequent [nrg++] = *tabname[0];}
  
  tmpstr.erase(tmpstr.begin(),tmpstr.end());
  tmpstr += *(tabname[0]);
  
  size_t st = tabname[0]->size();
  
  for (int j = 1; j < size; j++)
  {
    for (int l =1; l < size ; l++)
    {
      if (l != j) tmpstr += *(tabname[l]);     
    }
    double indic = conftwo(mappy[tmpstr], mappy[*tabname[j]],cursup);
    if (indic > Minconf)  
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
      allrules(antecedant[m],consequent[m],consalpha,nrg,m+1,mappy,cursup,size-1,Minconf, ruliste,strtemp);
    }
  }    
  
  if (tree.son) {genrules(*tree.son,size+1,tabname,strtemp,varnames,mappy,Minconf,ruliste);}
  if (tree.brother) {genrules(*tree.brother,size,tabname,str,varnames,mappy,Minconf,ruliste);}
  
  delete [] antecedant;
  delete [] consequent;
  delete [] consalpha; 
} 


void Genrulesroot (pnodesr & roots, int size, std::string ** tabname, std::string &str, std::string * varnames, std::unordered_map<std::string,double>& mappy,double & Minconf, std::list<rules>& ruliste)
{ std::string strtemp = str + varnames[roots.litem];
  double cursup = roots.sup;
  mappy[strtemp] = cursup;
  tabname[0] = &varnames[roots.litem];
  if (roots.son) genrules (*roots.son,size+1,tabname,strtemp,varnames,mappy, Minconf,ruliste);
  if (roots.brother) Genrulesroot(*roots.brother,size,tabname,str,varnames,mappy, Minconf,ruliste);    
}


int main (int argc , char ** argv)
{
  double * supvalue = NULL;
  int * sizevalue = NULL;
  int * litemvalue = NULL;
  std::vector<std::string> varnames;
  int nvar =0; int nrows = 0;
  char  * pathr = (char*)argv[1]; //  path to file of info freqset
  char  * pathv = (char*)argv[2];
  std::string confparam ((char*)argv[3]);
  
  if (confparam[0]=='c' && confparam[1]=='=')
      {
        confparam.erase(confparam.begin(),confparam.begin()+2);
      }
  else { std::cout << "Wrong third argument, please set c=yourMinconf" << std::endl;
          return 0;}

  double minConf = std::stod(confparam);

  std::string nbcoeffparam ((char*)argv[4]);

    if (nbcoeffparam[0]=='n' && nbcoeffparam[1]=='=')
      {
        nbcoeffparam.erase(nbcoeffparam.begin(),nbcoeffparam.begin()+2);
      }
  else { std::cout << "Wrong fourth argument, please set n=number_of_supplementary_coeff" << std::endl;
          return 0;}
  int nbcoeff = std::stoi(nbcoeffparam);
  
  int b = 5;
  std::string coeffstr [nbcoeff]; 

  for (int a = 0; a < nbcoeff;a++,b++)
    {
      coeffstr [a] = std::string ((char*)argv[b]);
    }


 init_data(pathr, pathv, nrows,nvar, varnames,supvalue,sizevalue,litemvalue);

  std::unordered_map<std::string,int> mapparam;
  mapparam["conftwo"] = 0;
  mapparam["cov"] = 1;
  mapparam["corr"] = 2;
  mapparam["kappa"] = 3;
  mapparam["mp"] = 4;

  std::cout << "Build PrefRecTRee_Rules ... "; 


  
  pnodesr * rootrules = new pnodesr(0,0,0);
  pnodesr* tabpnodes [nrows];
  gen_listetree(supvalue,litemvalue,sizevalue,nrows,tabpnodes);
  int sit =0;
  int k = nrows /10000;
  int restek = nrows -10000*k;
  int ct = 1;

  int nbg = 0;
  
  clist = tabpnodes[0];
  ctree = rootrules;
  cfath = rootrules;
  for (int j = 0; j < k; j++, ct++)
    { nbg = 10000*ct; 
      gen_tree (clist,ctree,cfath,tabpnodes,sit,nbg); 
    } 
  gen_tree(clist,ctree,cfath,tabpnodes,sit,nrows);

  std::cout << "done" << std::endl;

  std::unordered_map<std::string,double> mappysup;
  mappysup.reserve(2*nrows);
 
  std::list<rules> tabofrules;

  std::cout << "Testing " << "confidence" << " ... ";
  std::string ** tabname = new std::string * [varnames.size()];
  std::string st;
  st.erase(st.begin(),st.end());
  st.reserve (50);
  Genrulesroot (*rootrules->son,1,tabname,st,&varnames[0],mappysup, minConf ,tabofrules);
  std::cout << "done" << std::endl;
  std::cout << " There is " << tabofrules.size() << " confident rules " << std::endl;
  int nbrules = tabofrules.size();
  std::cout << "done" << std::endl;
   
  std::cout << "Calculate wished coefficient . . . " << std::endl;
  std::list<rules>::iterator ite = tabofrules.end();
  std::vector<std::string> antetab (nbrules);
  std::vector<std::string> comptab (nbrules);
  std::vector<std::string> settab (nbrules);
  std::vector<double> valuetab (nbrules);
  std::vector<std::vector<double>> tabcoeff (nbcoeff);
  for (int c = 0; c < nbcoeff; c++)
    {
      tabcoeff[c] = std::vector<double> (nbrules);
      double * ptabcoeff = &tabcoeff[c][0];
    }

  int i =0;
  for (std::list<rules>::iterator it = tabofrules.begin(); it!= ite;it++)
    {
      antetab[i] = it->ante;
      comptab[i] = it->comp;
      settab[i] = it->set;
      valuetab[i++] = it->conf;
    }

   fctionpt  dico;
 
  int par = 0;
  for (int j =0; j < nbcoeff;j++)
    {
      par = mapparam[coeffstr[j]];
      switch (par)
      {
        case 0:
        dico.mafonc = &conftwo; 
        break;
        case 1:
        dico.mafonc = &cov;
        break;
        case 2:
        dico.mafonc = &corr;
        break;
        case 3:
        dico.mafonc = &kappa;
        break;
        case 4:
        dico.mafonc = &MaxwellP;
        break;
      }

   double * ptabcoeff = &tabcoeff[j][0];
   for (int h = 0; h < nbrules; h++)
      {
        ptabcoeff[h] = dico.mafonc(mappysup[antetab[h]],mappysup[comptab[h]],mappysup[settab[h]]);
      }
 }
  std::cout << "done" << std::endl;
  free(supvalue); free(sizevalue); free (litemvalue);
  std::cout << "Write in file ..." << std::endl;

  if (argc > 5 + nbcoeff)
  {  
  char  * argpathout = (char*)argv[argc-1];
  std::string pathout (argpathout);
  pathout += ".txt";
  std::ofstream outfluxrules;
   outfluxrules.open(std::string(pathout));
   outfluxrules << "antecedant\t" << "complementary\t" << "confidence";
   int l = 0;
    for (l = 0; l < nbcoeff;l++)
      {
        outfluxrules << "\t" << coeffstr[l];
      } 
    outfluxrules << "\n"; 
    for(i = 0; i < nbrules ;i++)
      {
        outfluxrules << antetab[i] << "\t" << comptab[i] << "\t" << valuetab[i];
        for (l=0; l<nbcoeff;l++)
          { 
            outfluxrules << "\t" << tabcoeff[l][i] ;
          }
        outfluxrules << "\n";  
      }
    outfluxrules.close();
  }

  
  return 0;
}
 


