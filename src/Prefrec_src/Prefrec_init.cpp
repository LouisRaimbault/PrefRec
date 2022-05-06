#include "Prefrec_init.h"


const size_t SUL = sizeof(uint64_t);
const size_t SI = sizeof(int);
const size_t SPN = sizeof(pnodes);
const size_t SRN = sizeof(rnodes);


void Init_data (char * pathfile, int & nrows, int & maxul ,int & nvar, uint64_t **& Bitdata, std::vector<std::string> & varnames, char deli, double d_sup_init ) 
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
  std::cout << "nrow = " << nrows << "\n";
  for (t =0; t < nrows; t++, sit++ )
    { while (tabchar[sit] != '\n') {st += tabchar[sit++];}
      transaction[t]= st;
      st = "";
    }
  free (tabchar);
  int sup_init = nrows * d_sup_init;
  std::unordered_map<std::string,uint64_t*> maptr;
  std::unordered_map<std::string,int> map_int;
  st = "";
  int ko = 0;

  for (ko = 0; ko < nrows; ko++)
    { 
      for (char c :transaction[ko]) 
        {
          if (c!= deli) {st+= c;}            
          else {if(!map_int[st]) {map_int[st] = 0;}; map_int[st]++; st="";}
        }
         if (!map_int[st]) {map_int[st] = 0;} map_int[st]++; st="";
    }
  std::cout << "Done. Preparing BitData store in uint64_t ... \n";
  nvar = 0;
  maxul = nrows/64;
  int rst = nrows-64*maxul;
  if (rst) maxul++;
  std::unordered_map<std::string,int>::iterator it_int;
  std::unordered_map<std::string,int>::iterator it_itend = map_int.end();
  for (it_int = map_int.begin(); it_int != it_itend; it_int++)
    {
      if (it_int->second >= sup_init) {nvar++;}
    }

  Bitdata = (uint64_t**) malloc (nvar*sizeof(uint64_t*));
  std::cout << "nvar = " << nvar << "\n";
  if (nvar > 5000) {std::cout << "to much var \n"; }
  for (t = 0; t < nvar;t++) {Bitdata[t] = (uint64_t*)calloc(maxul,SUL);}

  varnames.reserve(nvar+1);
  std::unordered_map<std::string,uint64_t*>::iterator itmap;
  std::unordered_map<std::string,uint64_t*>::iterator itendmap = maptr.end();
  it_itend = map_int.end();
  t = 0;
  for (it_int = map_int.begin(); it_int != it_itend; it_int++)
    {
      if (it_int->second < sup_init) {continue;} 
      
      varnames.push_back(it_int->first+',');
      maptr[it_int->first] = Bitdata[t++]; 
    }

  t=0; i = 0; sit = 0;


  t = 0;
  for (ko = 0; ko < nrows;ko++)
  { st ="";
    for (char c :transaction[ko])
      { if (c!= deli) {st+= c;}
        else {if(map_int[st]>=sup_init) {maptr[st][t]|= 1UL << i;} st="";}
      }

    if(map_int[st]>=sup_init) {maptr[st][t]|= 1UL << i;}
    sit++; 
    t = sit/64;
    i = sit%64;
  }
  map_int.clear();
  maptr.clear();
  std::cout << " ok \n";
}



void sortoutdata (uint64_t** Bdata, int * sumvec, int nb_elem, std::vector<std::string>& varnames, char ordre)
{
  int i,j,treme,itreme,temp;
  uint64_t * tempt;
  std::string tempstr;
  if (ordre=='d')
  {  
    for (i = 0; i<nb_elem-1;i++)
    {itreme = i; treme=sumvec[i];
    for (j=i+1;j<nb_elem;j++) {if (sumvec[j]>treme) {treme = sumvec[j];itreme=j;}}
    temp=sumvec[itreme];sumvec[itreme]=sumvec[i];sumvec[i]=temp;
    tempt = Bdata[itreme]; Bdata[itreme] = Bdata[i];Bdata[i]=tempt;
    tempstr = varnames[itreme]; varnames[itreme] = varnames[i];varnames[i]=tempstr;
    }
  }
  if (ordre == 'i')
  {  
    for (i = 0; i<nb_elem-1;i++)
    {itreme = i; treme=sumvec[i];
    for (j=i+1;j<nb_elem;j++) {if (sumvec[j]<treme) {treme = sumvec[j];itreme=j;}}
    temp=sumvec[itreme];sumvec[itreme]=sumvec[i];sumvec[i]=temp;
    tempt = Bdata[itreme]; Bdata[itreme] = Bdata[i];Bdata[i]=tempt;
    tempstr = varnames[itreme]; varnames[itreme] = varnames[i];varnames[i]=tempstr;
    }
  }
}

int init_prefixtree (uint64_t ** Bitdata, uint64_t **& Bdata,std::vector<std::string>& varnames, int support, int nvar, int maxul, char ordre )
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
    if (ordre == 'i' || ordre == 'd') 
    { int sumvecord [tokeep];
      for (i = 0; i < nvar;i++) {if (sumvec[i] > support) {Bdata[j] = Bitdata[i]; sumvecord[j] = sumvec[i]; varnames[j++] = varnames[i];}} 
      sortoutdata(Bdata, sumvecord,tokeep,varnames,ordre);
      return tokeep;
    }
    for (i = 0; i < nvar; i++) {if (sumvec[i]>support) {Bdata[j] = Bitdata[i]; varnames[j++] = varnames[i];} }
    
  return tokeep;
}

void erase_freq (pnodes* Tree)
{ 
 if (Tree->son) {erase_freq(Tree->son);}
 if (Tree->brother) {erase_freq(Tree->brother);}
 free(Tree);

}


void extract_and_erase_freq (pnodes* Tree, std::string str, std::string *  listenom, double nrs , int size ,std::ofstream  & flux_set )
{ 
 
 std::string strtemp = str + listenom[Tree->litem];
 flux_set << strtemp << "\t" << Tree->sup << "\t" <<(double)Tree->sup/nrs << "\t" << size <<"\n"; 

 if (Tree->son) {extract_and_erase_freq(Tree->son,strtemp,listenom,nrs,size+1,flux_set);}
 if (Tree->brother) {extract_and_erase_freq(Tree->brother,str,listenom,nrs,size,flux_set);}
 free(Tree);
}


void extract_and_erase_freq (pnodes* Tree, std::string str, std::string *  listenom, double nrs, int size, std::ofstream  & flux_set, std::ofstream  & flux_rules )
{
 std::string strtemp = str + listenom[Tree->litem];
 double rel_sup = (double)Tree->sup/nrs;
 flux_set << strtemp << "\t" << Tree->sup << "\t" <<rel_sup << "\t" << size <<"\n"; 
 flux_rules << rel_sup << "\t" << size << "\t" << Tree->litem << "\n";  

 if (Tree->son) {extract_and_erase_freq(Tree->son,strtemp,listenom,nrs,size+1,flux_set,flux_rules);}
 if (Tree->brother) {extract_and_erase_freq(Tree->brother,str,listenom,nrs,size,flux_set,flux_rules);}
 free(Tree);
}