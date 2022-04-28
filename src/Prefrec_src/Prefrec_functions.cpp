#include "Prefrec_functions.h"


#pragma optimization_level 3
#pragma GCC optimize("Ofast,no-stack-protector,unroll-loops,fast-math,O3")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx")
#pragma GCC optimize("Ofast")
#pragma GCC target("avx,avx2,fma")
#pragma GCC optimization ("unroll-loops")

int nb_freq ;
int curlitem;
int minsup;
int nbc;
 static uint64_t Ultab [65];
 uint64_t * PT_ULTAB = &Ultab[0];




void creabitfield ( int *ind,int* mod,uint64_t * ptul, uint64_t * nlgtabl,  int nbca)
{   auto * yb = &Ultab[0]; uint64_t u; int t = 0; auto * ybi = &Ultab[0];
    auto * yd = &Ultab[64];
    for (auto a =0; a < nbca ; a++)
    {  
      for (yb = ybi; yb != yd; yb++)
      { 
        u ^= (-((ptul[ind[t]]>>mod[t++]) & 1UL) ^ u) & *yb;
      }
      nlgtabl[a] = u;
    }
}

void creabitfieldr ( int *ind,int * mod,uint64_t * ptul, uint64_t * nlgtabl,  int nbdases,  int &reste)
{   
    auto * yb = &Ultab[0]; uint64_t u; int t = 0;
    auto * yd = &Ultab[64]; auto * ybi = &Ultab[0];
    for (auto a =0; a < nbdases ; a++)
    {  
      for (yb = ybi; yb != yd; yb++)
      { 
        u ^= (-((ptul[ind[t]]>>mod[t++]) & 1UL) ^ u) & *yb;
      }
      nlgtabl[a] = u;
    }
    u = 0;
    for ( yb = ybi; yb != &Ultab[reste]; yb++) 
    {
     u ^= (-((ptul[ind[t]]>>mod[t++]) & 1UL) ^ u) & *yb;}
    nlgtabl[nbdases] = u;

}



void init_Rnodes (rnodes* outrnodes, int* vecsom,  int* tindices ,int * mod,uint64_t** Bitdata, int &nbleft, int& nbcases, int Mins,  int &reste)
{ 
  
  int sum = 0; 
  if (reste ==0)
    {for (auto i=0; i < nbleft;i++)
      { sum = vecsom[i];
        if (sum<minsup) {outrnodes[i].tab=NULL; outrnodes[i].indic = 0;continue;}
        nb_freq++;
        uint64_t * ptu = (uint64_t*) malloc(nbcases*SUL);
        creabitfield(tindices,mod,Bitdata[i],ptu,nbcases);
          outrnodes[i].sup = sum;
          outrnodes[i].tab=ptu;
          outrnodes[i].indic = 1;
      }
    }

  else {
        for (auto i=0; i < nbleft;i++)
          { sum = vecsom[i];
            if (sum<minsup) {outrnodes[i].tab=NULL; outrnodes[i].indic = 0;continue;}
            nb_freq++;
            uint64_t * ptu = (uint64_t*) malloc(nbcases*SUL);
            creabitfieldr(tindices,mod,Bitdata[i],ptu,nbcases-1,reste);
            outrnodes[i].sup = sum;
            outrnodes[i].tab=ptu;
            outrnodes[i].indic = 1;
          }
      }
}



void depthwalk (pnodes * Tree, uint64_t * cand, rnodes* tabrnodes)
{ int elem = Tree->litem;
  bool d = !Tree->brother;
  if (!tabrnodes[elem].indic) {if (!d) {depthwalk(Tree->brother, cand, tabrnodes);} return;}
  int supcand = 0;
  auto * tib = tabrnodes[elem].tab;
  auto * ed = &cand[nbc];
  for (auto * tob = cand; tob != ed; ++tob,++tib )
    {supcand += __builtin_popcountl(*tob&*tib);}
  if (supcand<minsup) 
  {tabrnodes[elem].indic=0;
  if (!d) depthwalk(Tree->brother, cand, tabrnodes);
  tabrnodes[elem].indic=1;
  return;
  }
           
  nb_freq++;
  if (!Tree->son)
    {Tree->son = (pnodes*)malloc(SPN);
     Tree->son->sup = supcand;
     Tree->son->litem = curlitem;
     Tree->son->son = NULL;
     Tree->son->brother = NULL;
     if (!d) {depthwalk(Tree->brother, cand, tabrnodes);}
     return;
    }   
  pnodes * trpnodes = Tree->son;
  uint64_t  * ptab = (uint64_t*) malloc (nbc*sizeof(uint64_t)) ;
  auto * pptab = &ptab[0];
  tib = tabrnodes[elem].tab;
  for (auto * tob = cand; tob != ed; ++tob,++tib,++pptab ){*pptab = *tob&*tib;}        
  depthwalk(trpnodes,ptab,tabrnodes);
  free(ptab);
  Tree->son = (pnodes*)malloc(SPN);
  Tree->son->sup = supcand;
  Tree->son->litem = curlitem;
  Tree->son->brother = trpnodes;
  Tree->son->son = NULL; 
  if (!d) depthwalk(Tree->brother, cand, tabrnodes);
  return;       
}


void rootwalk (pnodes * root ,rnodes * rootnodes,  int nbcases, int nbleft )
{ pnodes * temp = root->brother;
  pnodes * trpn = NULL;
  for (int i =nbleft-1; i > -1 ; i--)
    {if (rootnodes[i].indic)
        {   trpn = temp->son;
            if ( trpn == NULL) 
              {temp->son = (pnodes*)malloc(SPN);
               temp->son->sup = rootnodes[i].sup,temp->son->litem = curlitem; temp->son->brother = NULL; temp->son->son = NULL;
              } 
            else 
            { 
              depthwalk(trpn,rootnodes[i].tab,rootnodes);
              temp->son = (pnodes*)malloc(SPN);
              temp->son->sup = rootnodes[i].sup,temp->son->litem = curlitem; temp->son->brother = trpn; temp->son->son = NULL;
            }
            
        }
      temp = temp->brother;
    }
}



void Bitprefrec (uint64_t** Bitdata, std::vector<std::string>& varnames, int maxul, int nvar, char ordre, pnodes *& root)
{  curlitem = 0; nb_freq = 0;
 std::cout << "Support value is " << minsup << ". Checking unfrequent items ... \n ";
 int sum = 0;
 auto a = 0; int b= 0;
 auto e = 0;
 auto n = 0;
 uint64_t ** Bdata = NULL ;
 auto * pul = &Bitdata[0][0];
 auto * pulb = &Bitdata[0][0];
 auto * pulc = &Bitdata[0][0];
 int keep = init_prefixtree(Bitdata, Bdata ,varnames, minsup,nvar,maxul,ordre);
 std::cout <<"Done The Supportvalue is " << minsup << ". Start PrefRec with  " << keep << " frequents variables ..." <<  std::endl;
 if(keep ==0)
  {
    std::cout << "Not a single frequent itemSet, please try with a lower relative minSup value. " << std::endl;
  }
 
 auto start = std::chrono::system_clock::now();
 root = (pnodes*)malloc(SPN);
 root->son = NULL;
 root->brother = NULL;
 root->litem = 0;
 root->sup = 0;
 for (int j=0; j<keep;j++)
  { 
   sum = 0;
   uint64_t u;
   curlitem = j;
   nb_freq++;
   pul = Bdata[j];
   pulc = &Bdata[j][maxul];
   for (a = 0; a < maxul;a++)
   {sum+= __builtin_popcountl(pul[a]);}

   nbc = sum/64;
   int reste = sum-64*nbc;
   if (reste) nbc ++;
   e = 0;
   rnodes tabrn [j+1]; 
   int * indices = (int*) malloc(sum*sizeof(int));
   int * mod = (int*) malloc(sum*sizeof(int));
   for (a=0; a < maxul; a++)
    { u = pul[a];
      for (n = 0; n < 64; n++)
       { if ((u >> n) & 1UL) {indices[e]= a; mod[e++]=n;} }
    }

    int rootsum [j+1]{0};
    int s = 0;
    for (a = 0; a < j; a++)
      {
        s = 0; pulb = Bdata[a];
        for (pul = Bdata[j]; pul != pulc; ++pul,++pulb)
          { s += __builtin_popcountl(*pul & *pulb);}
        rootsum[a] = s;
      }

  init_Rnodes(tabrn,rootsum,indices,mod,Bdata,j,nbc,minsup,reste);

  rootwalk(root,tabrn,nbc,j);

  pnodes * nouvnod = (pnodes*)malloc(SPN);
  nouvnod->sup = sum;
  nouvnod->litem = curlitem;
  nouvnod->brother = root->brother;
  nouvnod->son = NULL;
  root->brother = nouvnod;
  for (a = 0; a < j; a++)
   {if (tabrn[a].tab) free(tabrn[a].tab);}
   free(indices);
   free(mod);
  }

  std::chrono::duration<double> diffe= std::chrono::system_clock::now() - start;
  std::cout << "End PrefRec, there is  " << nb_freq <<" sets, extracted int  " <<diffe.count() << " secs"<< std::endl;
  free (Bdata);

}