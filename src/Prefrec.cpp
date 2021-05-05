#include "Prefrec.h"
#pragma optimization_level 3
#pragma GCC optimize("Ofast,no-stack-protector,unroll-loops,fast-math,O3")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx")
#pragma GCC optimize("Ofast")//Comment optimisations for interactive problems (use endl)
#pragma GCC target("avx,avx2,fma")
#pragma GCC optimization ("unroll-loops")


dimentions::dimentions(int nb_elem_, int nb_lignes_) : nb_elem (nb_elem_) , nb_lignes (nb_lignes_) {};

dimentions longueurfichier (char * cheminfichier) //détermine la longueur d'un fichier 
{ 
  FILE *Fichier = fopen(cheminfichier,"r"); // ouvre le fichier en mode lecture 
  int longueur = 0;
  int nb_lignes = 0;
  char c;

  if (Fichier != NULL)
    { 
      while (fscanf(Fichier,"%c",&c)!= EOF)  // Tant que le pointeur dans le fichier n'est pas sur l'end of file
       {longueur ++;
        if (c == '\n') nb_lignes++; }                      

    }
    fclose(Fichier);
    dimentions out (longueur, nb_lignes); 
    return out;  

}

void remplirfichier (char *cheminfichier, char* chaine, dimentions& dimtableau) // prends une chaine de character ayant la bonne chaine 
{ 
  FILE * Fichier = fopen(cheminfichier,"r");
  if ( Fichier != NULL && chaine  != NULL)
    { int i = 0;
      char c;
     while (fscanf(Fichier , "%c",&c)!= EOF) // tant que exctraction de la donnée character  dans la valeure c différente de end of file
      {chaine[i]=c;
        i++;}}
        else perror ("\n\n memorisierfichier ");

      fclose(Fichier);     
}

void transformintotransac (char* chainecara, std::string* outtransac,dimentions& dimtab)
{ 
  int nbligne = dimtab.nb_lignes;
  int sit = 0;
  std::string aremp ="";
  aremp.reserve(100);
  for (int i =0; i < nbligne; i++, sit++ )
   { while (chainecara[sit] != '\n')
     {aremp += chainecara[sit++];}
      outtransac[i]= aremp;
      aremp = "";}

}

void getmap  (std::vector<std::string> transac, std::unordered_map<std::string,int>& maptr,std::vector<std::string> &listename ,char deli)
{ 
  
  std::string qs;
  qs.reserve(100);
  for (std::string s :  transac)
  { 
    for (char c :s)
      { if (c!= deli) {qs+= c;}
      else {maptr[qs]; qs.erase(qs.begin(),qs.end());}
      }
    maptr[qs]; qs.erase(qs.begin(),qs.end());
  }
  
}

void transactiontoBitmax( bool ** mytab,std::unordered_map<std::string,int> &mappy,std::vector<std::string> transac, std::vector<std::string> &listename ,char deli,int nb_individus)
{ 
  std::unordered_map<std::string,int>::iterator itmap;
  std::unordered_map<std::string,int>::iterator itendmap = mappy.end();
  int a = 0;
  for (itmap = mappy.begin();itmap != itendmap;itmap++)
  {itmap->second = a++;
    listename.push_back(itmap->first);}
  int i = 0;
  for (std::string s :  transac)
  { std::string q ="";
    for (char c :s)
      { if (c!= deli) {q+= c;}
      else {mytab[mappy[q]][i]=1; q.erase(q.begin(),q.end());}
      }
    mytab[mappy[q]][i++]=1;
  }
  
  
}



pnodes::pnodes (int sup_, unsigned int litem_, pnodes* brother_ ):sup(sup_), litem(litem_), brother(brother_) {};


void tri_tableau (bool** dataframept, std::vector<int> mytab, int nb_elem, std::vector<std::string> &listenoms)
{
  int i,j,min,imin,temp;
  bool* tempt;
  std::string tempstr;
  for (i = 0; i<nb_elem-1;i++)
  {imin = i; min=mytab[i];
  for (j=i+1;j<nb_elem;j++) {if (mytab[j]>min) {min = mytab[j];imin=j;}}
  temp=mytab[imin];mytab[imin]=mytab[i];mytab[i]=temp;
  tempt = dataframept[imin]; dataframept[imin] = dataframept[i];dataframept[i]=tempt;
  tempstr = listenoms[imin]; listenoms[imin] = listenoms[i];listenoms[i]=tempstr;
  }
}

int init_prefixtree (bool** mydataframe,std::vector<std::string>& listenoms, int support, int nbvar, int nbind )
{ 
  std::vector<int> somme_vec (nbvar);
  
  int tokeep = 0;
  int sumb = 0;
  bool * botab = NULL;
  int i = 0;

  for (int h = 0; h <nbvar; h++ )
  { botab = mydataframe[h];
    sumb = 0;
    for (i = 0; i < nbind;i++)
      {
        sumb += botab[i];
      }
    somme_vec[h]=sumb;
    if (sumb > support ){tokeep++;}}
  tri_tableau(mydataframe, somme_vec,nbvar,listenoms);
  return tokeep;
}

void Bitmaxtranspos (bool ** Bitranspose,bool ** Bitmax, int nbc, int nbl)
{ 
  for (int i=0; i < nbc; i++)
  {Bitranspose[i]= new bool [nbl];
    bool * addelem = Bitranspose[i];
    for (int j=0; j < nbl; j++)
    {addelem[j]=Bitmax[j][i];}
    }
}


rnodes::rnodes ( int sup_,bool indic_, uint64_t* tab_):sup(sup_), indic(indic_), tab(tab_){};


int nbfreq = 0;
double minsup =0;
unsigned int curlitem = 0;


void rootsumtest ( int * outsum, int* indyces, bool ** TRBitmax , int supcand, int nbleftrn)
{
    for (int i=0; i < supcand; i++)
        {bool * ptshT = TRBitmax[indyces[i]];
          for (int j=0; j < nbleftrn;j++)
           {outsum[j] += ptshT[j];}
        }     
}

void transformintolong ( int *ind,bool * ptshort, uint64_t * nlgtabl,  int nbc)
{   int y =0;
    std::bitset<64> bset;
    for (int a =0; a < nbc ; a++)
    {
      for (y =0; y < 64; y++)
      {bset[y]=ptshort[*ind++];}
      nlgtabl[a] = bset.to_ulong();}
}

void transformintolongreste ( int *ind,bool * ptshort, uint64_t * nlgtabl,  int nbdases,  int &reste)
{   int y =0;
    std::bitset<64> bset;
    for (int a =0; a < nbdases ; a++)
    {
      for (y =0; y < 64; y++)
      {bset[y]=ptshort[*ind++];}
      nlgtabl[a] = bset.to_ulong();}
    std::bitset<64> asetone;
    for ( y = 0; y < reste; y++) {asetone[y]= ptshort[*ind++];}
    nlgtabl[nbdases] = asetone.to_ulong();

}


void createvecrnodes (rnodes* outrnodes, int* vecsom,  int* tindices ,bool** Bitmax, int &nbleft, int& nbcases, int Mins,  int &reste)
{ 
  int sit = 0;
  int sum = 0; 
  if (reste ==0)
    {for (int i=0; i < nbleft;i++)
      { sum = vecsom[i];
        if (sum>Mins)
          { nbfreq++;
            uint64_t * ptu = new uint64_t[nbcases];
            transformintolong(tindices,Bitmax[sit],ptu,nbcases);
            outrnodes[i]= rnodes(sum,1,ptu);
          }
        else outrnodes[i]=rnodes(0,0,NULL);  
           sit++;
      }
    }

  else {
        for (int i=0; i < nbleft;i++)
          { sum = vecsom[i];
            if (sum>Mins)
            { nbfreq++;
              uint64_t *ptu  = new uint64_t [nbcases];
              transformintolongreste(tindices,Bitmax[sit],ptu,nbcases-1,reste);
              outrnodes[i]= rnodes(sum,1,ptu);
            }
            else outrnodes[i]=rnodes(0,0,NULL);      
           sit++;
          }
      }
}

void depthwalk (pnodes & Tree, uint64_t * cand, rnodes* tabrnodes , int nbcases)
{ unsigned int elem = Tree.litem;
    if (tabrnodes[elem].indic)
      {   
        int supcand = 0;
        uint64_t * tib = tabrnodes[elem].tab;
        for (int i =0; i < nbcases;i++)
          {supcand += __builtin_popcountl(cand[i]&tib[i]);}
        if (supcand>minsup)
          {   
            nbfreq++;
            if (Tree.son == NULL) Tree.son = new pnodes(supcand, curlitem,NULL);
            else 
                {
                  pnodes * trpnodes = Tree.son;    
                  uint64_t * ptab = new uint64_t [nbcases];
                  uint64_t * ntab = tabrnodes[elem].tab;
                  for (int i =0; i < nbcases;i++)
                    {
                      ptab[i] = cand[i]&ntab[i];
                    }
                  depthwalk(*trpnodes,ptab,tabrnodes,nbcases);
                  Tree.son = new pnodes (supcand,curlitem,trpnodes);
                  delete [] ptab;
                }


            if (Tree.brother) depthwalk(*Tree.brother, cand, tabrnodes,nbcases);
        }
                
        else 
            { if (Tree.brother)
               {  
                   tabrnodes[elem].indic=0;
                   depthwalk(*Tree.brother, cand, tabrnodes,nbcases);
                   tabrnodes[elem].indic = 1;
               }
            }     
        }

      else { if (Tree.brother)depthwalk(*Tree.brother, cand, tabrnodes,nbcases);}
}


void rootwalk (pnodes & root ,rnodes * rootnodes,  int nbcases, int nbleft )
{ pnodes * temp = root.brother;
  pnodes * trpn = NULL;
  for (int i =nbleft-1; i > -1 ; i--)
    {if (rootnodes[i].indic)
        {   trpn = temp->son;
            if ( trpn == NULL) temp->son = new pnodes (rootnodes[i].sup,curlitem,NULL);
            else 
            {  
              depthwalk(*trpn,rootnodes[i].tab,rootnodes,nbcases);
              temp->son = new pnodes (rootnodes[i].sup,curlitem,trpn);
            }
            
        }
      temp = temp->brother;
    }
}


pnodes Bitprefrec (bool** Bitmax, std::vector<std::string> &varnames, double relativeSup, int nbind, int nbvar)
{
 minsup = nbind* relativeSup -1;
 nbfreq=0;
 curlitem = 0;
 int keep = init_prefixtree(Bitmax, varnames, minsup,nbvar,nbind);
  std::cout <<"The Supportvalue is " << minsup+1  << std::endl;
  if(keep ==0) {
    std::vector<std::string> outempty;
    outempty.push_back("empty");
    std::cout << "Not a single frequent itemSet, please try with a lower relative minSup value. " << std::endl;
    return pnodes (0,0,NULL);
  }
 
    bool * TBitmax [nbind]; 
  Bitmaxtranspos(TBitmax,Bitmax,nbind,keep);
  std::cout << "Start PrefRec with  " << keep << " frequents variables ... "<<  std::endl;
   clock_t stimer, etimer;
   stimer = clock();
   pnodes rootalpha (0,0,NULL);

   for (int j=0; j<keep;j++)
 { 
   curlitem = j;
   nbfreq++;    
   bool * binarynewvec = Bitmax[j];
   int sum = 0;
   int a = 0;
   for (a = 0; a < nbind;a++)
   {sum+=binarynewvec[a];}

  
   int* indices = new int[sum];
   int e = 0;
   rnodes * tabrn =  new rnodes [j+1]; 
   for (int a=0; a < nbind; a++)
    {if (binarynewvec[a]) indices[e++]=a; }

    int nbc = sum/64;
    int reste = sum-64*nbc;
    int* rootsom = new int [j+1]{0};

    rootsumtest(rootsom,indices,TBitmax,sum,j);

   if (reste != 0) nbc ++;
  createvecrnodes(tabrn,rootsom,indices,Bitmax,j,nbc,minsup,reste);
  rootwalk(rootalpha,tabrn,nbc,j);

  pnodes * nouvnod = new pnodes (sum,curlitem,NULL);
  nouvnod->brother = rootalpha.brother;
  rootalpha.brother = nouvnod;


  for (a = 0; a < j; a++)
   {if (tabrn[a].tab) delete [] tabrn[a].tab;}
  delete [] indices;
  delete [] rootsom;
  delete [] tabrn;
   }
 etimer = clock();
double diffe = double(etimer - stimer)/double(CLOCKS_PER_SEC);
 std::cout << "End PrefRec, there is  " << nbfreq <<" sets, extracted int  " <<diffe << " secs"<< std::endl;

 return rootalpha;
    
}

void extract_and_erase_freq (pnodes*& Tree, std::string str, std::string *  listenom, std::string * namevalue, int * supvalue, int* sizevalue, int* litemvalue , int size ,int& situeur )
{ 
 std::string strtemp = str + listenom[Tree->litem];
 namevalue[situeur] = strtemp;
 litemvalue[situeur] = Tree->litem;
 sizevalue [situeur] = size;
 supvalue[situeur++] = Tree->sup;
 
 if (Tree->son) {extract_and_erase_freq(Tree->son,strtemp,listenom,namevalue,supvalue,sizevalue,litemvalue,size+1,situeur);}
 if (Tree->brother) {extract_and_erase_freq(Tree->brother,str,listenom,namevalue,supvalue,sizevalue,litemvalue,size,situeur);}
delete Tree;

}

int main (int argc , char ** argv )
{

  std::cout << "Importing data from txt file . . . " ;  
  char  * monchemin = (char*)argv[1];
  dimentions dimtab = longueurfichier(monchemin);
  int nb_elem = dimtab.nb_elem;
  std::string delimparam ((char*)argv[2]);
  std::string relsupparam ((char*)argv[3]);
  char delim ;

  if (delimparam[0]=='d' && delimparam[1]=='=')
    {
      delim = delimparam[2];
    }
  else { std::cout << "Wrong second argument, please set d=yourdelim for separate item transactions" << std::endl;
         return 0;}

  if (relsupparam[0]=='s' && relsupparam[1]=='=')
    {relsupparam.erase(relsupparam.begin(),relsupparam.begin()+2);}
  else {std::cout << "Wrong third argument, please set s=yourrelativeminSup" << std::endl;
        return 0;}
  
  char *tabchar = new char [nb_elem];
  remplirfichier (monchemin,tabchar,dimtab);

  std::vector<std::string> vectransac (dimtab.nb_lignes);
  transformintotransac(tabchar,&vectransac[0],dimtab);
  
  int nb_individus = vectransac.size();
  std::vector<std::string> varnames;
  std::unordered_map<std::string,int> mappi ; 
  getmap (vectransac,mappi,varnames,delim);
  int nb_var = mappi.size();
  bool * Bitmax [nb_var];
  for (int s = 0; s < nb_var;s++)
    {Bitmax[s]= new bool [nb_individus]{0};} 
  transactiontoBitmax(Bitmax,mappi,vectransac,varnames,delim,nb_individus);   
  dimtab.nb_var = varnames.size();
  for (std::string & vst : varnames )
    {
      vst.insert(0,1,' ');
    }
  
   double relativeminsup = std::stod(relsupparam);
   std::cout << "ok " << std::endl;
   pnodes outnodes = Bitprefrec(Bitmax, varnames,relativeminsup,nb_individus,nb_var);

  if (nbfreq == 0) {return 0;}
  std::vector<int> Supvalue (nbfreq);
  std::vector<int> sizevalue (nbfreq);
  std::vector<int> litemvalue  (nbfreq);
  std::vector<double> RelativeSupvalue  (nbfreq);
  std::vector<std::string> namevalue (nbfreq);
  std::string ststring;
  ststring.reserve (50);
  int situ =0;
  double nbf = (double)nbfreq;
  extract_and_erase_freq (outnodes.brother, ststring, &varnames[0],&namevalue [0], &Supvalue[0], &sizevalue[0], &litemvalue[0],1,situ); 
  
  for (int o = 0; o < nbfreq; o++)
    {
      RelativeSupvalue[o] = double(Supvalue[o])/double(nb_individus);
    }

  std::cout << "writting in file ..." << std::endl;
  if (argc > 4)
    {
      char  * argoutpath = (char*)argv[4];
      std::string outpath (argoutpath);
      outpath += ".txt";
      std::ofstream outflux;
      outflux.open(std::string(outpath));
      outflux << "nameset\t" << "Support\t" << "relative_support"<<std::endl;
      for (int i = 0; i < nbfreq;i++)
        {
          outflux << namevalue[i] << "\t" << Supvalue[i] << "\t" <<RelativeSupvalue[i] << std::endl;
        }
      outflux.close();       
    }

  if (argc > 5)
    { int j=0;
      char  * argoutpathrules = (char*)argv[5];
      std::string outpathrules (argoutpathrules);
      outpathrules += ".txt";
      std::ofstream outfluxrules;
      outfluxrules.open(std::string(outpathrules));
      for(j = 0; j < nbfreq;j++)
        {
          outfluxrules << RelativeSupvalue[j] << "\t" << sizevalue[j] << "\t" << litemvalue[j] << std::endl; 
        }
      outfluxrules.close();
      std::string outpathvar (argoutpathrules);
      outpathvar += "_variables.txt";

      std::ofstream outvarflux;
      outvarflux.open(outpathvar);
      for(j = 0; j < nb_var;j++)
        {
          outvarflux << varnames[j] << std::endl;
        }   
      outvarflux.close();
    }

  std::cout << "done" << std::endl;
    return 0;
}
