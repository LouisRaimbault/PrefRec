#include "Prefrec_init.h"


const size_t SUL = sizeof(uint64_t);
const size_t SI = sizeof(int);
const size_t SPN = sizeof(pnodes);
const size_t SRN = sizeof(rnodes);


void Init_data (char * pathfile, int & nrows, int & maxul ,int & nvar, uint64_t **& Bitdata, std::vector<std::string> & varnames, char deli ) 
{ 
  std::cout << "import data from file as char array ... \n";
  char * tabchar;
  long longueur = 0;
  long i = 0;
  char c;
  int sit =0;
  int t = 0;
  std::string st ="";
  FILE *Fichier = fopen(pathfile,"r"); 

  if (Fichier != NULL)
    { 
      while (!feof(Fichier))  
       {c = getc(Fichier);
        longueur ++;
        if (c == '\n') nrows++; }                      

    }
  fclose(Fichier);
  tabchar = (char*)malloc(longueur * sizeof(char));
  Fichier = fopen(pathfile,"r");
  if ( Fichier != NULL && tabchar  != NULL)
    { 
     while (i < longueur)
      { c = getc(Fichier);
        tabchar[i++]=c;}
      }
  else perror ("\n\n problem in file ");
  fclose(Fichier);
  std::cout << " Done. Constructing as transasction ... \n";
  std::vector<std::string> transaction (nrows);      
  st.reserve(1000);
  for (t =0; t < nrows; t++, sit++ )
    { while (tabchar[sit] != '\n') {st += tabchar[sit++];}
      transaction[t]= st;
      st = "";
    }
  free (tabchar);

  std::unordered_map<std::string,uint64_t*> maptr;
  st = "";
  for (std::string s :  transaction)
    { 
      for (char c :s) 
        {
          if (c!= deli) {st+= c;}            
          else {maptr[st]; st="";}
        }
        maptr[st]; st="";
    }
  std::cout << "Done. Preparing BitData store in uint64_t ... \n";
  nvar = maptr.size();
  maxul = nrows/64;
  int rst = nrows-64/maxul;
  if (rst) maxul++;
  Bitdata = (uint64_t**) malloc (nvar*sizeof(uint64_t*));
  for (t = 0; t < nvar;t++) {Bitdata[t] = (uint64_t*)calloc(maxul,SUL);}

  varnames.reserve(nvar+1);
  std::unordered_map<std::string,uint64_t*>::iterator itmap;
  std::unordered_map<std::string,uint64_t*>::iterator itendmap = maptr.end();
  t=0; i = 0; sit = 0;
  for (itmap = maptr.begin();itmap != itendmap;itmap++)
  {
   itmap->second = Bitdata[t++];
   varnames.push_back(itmap->first+','); 
  }

  t = 0;
  for (std::string s : transaction)
  { st ="";
    for (char c :s)
      { if (c!= deli) {st+= c;}
        else {maptr[st][t]|= 1UL << i; st="";}
      }
    maptr[st][t]|= 1UL << i;
    sit++; 
    t = sit/64;
    i = sit%64;
  }

  std::cout << " ok \n";
}



void sortoutdata (uint64_t** Bdata, int * sumvec, int nb_elem, std::vector<std::string>& varnames, char order)
{
  
  if (order == 'a')
   { 
    int nnames = varnames.size();
    std::cout << "il y a " << nnames << "variables \n";
    int i,j,treme,itreme,temp;
    i = 0;
    uint64_t * tempt;
    std::string tempstr;
    int * corres_var = (int*)malloc(nnames*sizeof(int));

    for (std::string tata : varnames)
    { 
      tata.erase(tata.begin(),tata.begin()+1);
      tata.erase(tata.begin()+tata.size()-1); 
      corres_var[i] = std::stoi(tata);
      i++;   
    }

    for (i = 0; i<nb_elem-1;i++)
    {itreme = i; treme=corres_var[i];
    for (j=i+1;j<nb_elem;j++) {if (corres_var[j]<treme) {treme = corres_var[j];itreme=j;}}
    temp=corres_var[itreme];corres_var[itreme]=corres_var[i];corres_var[i]=temp;
    tempt = Bdata[itreme]; Bdata[itreme] = Bdata[i];Bdata[i]=tempt;
    tempstr = varnames[itreme]; varnames[itreme] = varnames[i];varnames[i]=tempstr;
    }

    free(corres_var);
  }
  
}

int init_prefixtree (uint64_t ** Bitdata, uint64_t **& Bdata,std::vector<std::string>& varnames, int support, int nvar, int maxul, char order )
{ 
  int sumvec [nvar];
  int tokeep = 0;
  int sumb = 0;
  uint64_t * ptbul = NULL;
  int i = 0; int j = 0;
  for (int h = 0; h <nvar; h++ )
  { ptbul = Bitdata[h];
    sumb = 0;
    for (i = 0; i < maxul;i++)
      {
        sumb += __builtin_popcountl(ptbul[i]);
      }
    sumvec[h]=sumb;
    if (sumb > support ) tokeep++;
  }
    Bdata = (uint64_t**) malloc(nvar*sizeof(uint64_t*));
       int sumvecord [tokeep];
      for (i = 0; i < nvar;i++) {if (sumvec[i] > support) {Bdata[j] = Bitdata[i]; sumvecord[j] = sumvec[i]; varnames[j++] = varnames[i];}} 
      sortoutdata(Bdata, sumvecord,tokeep,varnames,order);
      return tokeep;
    
    for (i = 0; i < nvar; i++) {if (sumvec[i]>support) {Bdata[j] = Bitdata[i]; varnames[j++] = varnames[i];} }
    
  return tokeep;
}

void extract_and_erase_freq (pnodes* Tree, int & nt)
{ nt++;
 if (Tree->son) {extract_and_erase_freq(Tree->son,nt);}
 if (Tree->brother) {extract_and_erase_freq(Tree->brother,nt);}
 free(Tree);
 Tree = NULL;

}


void erase_freq (pnodes* Tree)
{ 
 if (Tree->son) {erase_freq(Tree->son);}
 if (Tree->brother) {erase_freq(Tree->brother);}
 free(Tree);

}


void extract_and_erase_freq (pnodes* Tree, std::string str, std::string *  listenom, double nrs , int size ,std::ofstream  & flux_set )
{ 
 
 std::string strtemp = str + *Tree->link_rn->var;
 flux_set << strtemp << "\t" << Tree->sup << "\t" <<(double)Tree->sup/nrs << "\t" << size <<"\n"; 

 if (Tree->son) {extract_and_erase_freq(Tree->son,strtemp,listenom,nrs,size+1,flux_set);}
 if (Tree->brother) {extract_and_erase_freq(Tree->brother,str,listenom,nrs,size,flux_set);}
 free(Tree);
}


void extract_and_erase_freq (pnodes* Tree, std::string str, std::string *  listenom, double nrs, int size, std::ofstream  & flux_set, std::ofstream  & flux_rules )
{
 std::string strtemp = str + *Tree->link_rn->var;
 double rel_sup = (double)Tree->sup/nrs;
 flux_set << strtemp << "\t" << Tree->sup << "\t" <<rel_sup << "\t" << size <<"\n"; 
 flux_rules << rel_sup << "\t" << size << "\t" << Tree->link_rn->lab << "\n";  

 if (Tree->son) {extract_and_erase_freq(Tree->son,strtemp,listenom,nrs,size+1,flux_set,flux_rules);}
 if (Tree->brother) {extract_and_erase_freq(Tree->brother,str,listenom,nrs,size,flux_set,flux_rules);}
 free(Tree);
}




void set_rnodes_in_order (rnodes * root_rnodes, rnodes * n_rnodes)
{
  int cur_sup = n_rnodes->sup_data;
  while (root_rnodes->brother->sup_data < cur_sup ) { root_rnodes = root_rnodes->brother; }
  rnodes *  temp = root_rnodes->brother;
  root_rnodes->brother = n_rnodes;
  n_rnodes->brother = temp;
}

void erase_mobile (pnodes * pn, int * nbdel)
{
  if (pn->son) {erase_mobile(pn->son,nbdel);}
  if (pn->brother) {erase_mobile(pn->brother,nbdel);}
  (*nbdel)++;
  free (pn);
}


void erase_condition_father (pnodes * pn, rnodes  * n_rnodes, int * nbdel)
{ 
  pnodes * tpe_son = pn->son;
  pnodes * rt = pn->son;
  pnodes * temp = NULL;

  while (rt->brother)
    {
      if (rt->brother->link_rn == n_rnodes)
        {
          (*nbdel)++;
          if (rt->brother->son) {erase_mobile(rt->brother->son,nbdel);}
          temp = rt->brother;
          rt->brother = temp->brother;
          free(temp);
          rt = rt->brother;
          while (rt)
            {
              if (rt->son) {erase_condition_father(rt,n_rnodes,nbdel);}
              rt = rt->brother;
            }
          break;
        } 

      if (rt->brother->son) {erase_condition_father(rt->brother,n_rnodes,nbdel);}
      rt = rt->brother;     
    }

    if (tpe_son->link_rn == n_rnodes)
    {
      (*nbdel)++;
       pn->son = tpe_son->brother; 
       if (tpe_son->son) {erase_mobile(tpe_son,nbdel);}
       free(tpe_son);
    }  
}


void erase_var (rnodes * rn, int * nbdel)
{ 
  (*nbdel)++; 

  pnodes * del_tree = rn->cor_bsb->brother;

  pnodes * temp = del_tree->brother;
  
  if (temp) 
    {
      temp->link_rn->cor_bsb = rn->cor_bsb;
    }
  
  
  pnodes * ttemps = NULL;
  pnodes * rt = NULL;
  
  if (del_tree->son)
    {
      erase_mobile(del_tree->son,nbdel);
    }

  rn->cor_bsb->brother = temp;
  free(del_tree);

  while (temp)
    {
      if (temp->son) {erase_condition_father(temp,rn,nbdel);}
      temp = temp->brother;

    }

}

void erase_nouvar (pnodes * Tree, rnodes * rn ,int * nbdel)
{
  pnodes * temp = Tree->brother;
  while (temp)
    {
       if (temp->son) {erase_condition_father(temp,rn,nbdel);}
       temp = temp->brother;     
    }

}