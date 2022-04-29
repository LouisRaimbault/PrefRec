#include "Prefrec_functions.h"


#pragma optimization_level 3
#pragma GCC optimize("Ofast,no-stack-protector,unroll-loops,fast-math,O3")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx")
#pragma GCC optimize("Ofast")
#pragma GCC target("avx,avx2,fma")
#pragma GCC optimization ("unroll-loops")

int nb_freq ;
rnodes * curlitem;
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



void init_Rnodes (pnodes * Tree , int * ngroot , int* tindices ,int * mod, int& nbcases,  int &reste)
{ 
  int sit = 0;
  int sum = 0;
  Tree = Tree->brother; 
  if (reste ==0)
    { 
      while (Tree)
        { sum = Tree->link_rn->sup;         
          if (sum<minsup) {Tree->link_rn->sup = 0; Tree = Tree->brother; continue;}
          (*ngroot)++;
          uint64_t * ptu = (uint64_t*) malloc(nbcases*SUL);
          creabitfield(tindices,mod,Tree->link_rn->tab_data,ptu,nbcases);
          Tree->link_rn->sup = sum;
          Tree->link_rn->tab = ptu;
          Tree = Tree->brother;
        }
     nb_freq += (*ngroot);   
    }

  else {
        
        while (Tree)
          { sum = Tree->link_rn->sup;
            if (sum<minsup) {Tree->link_rn->sup = 0;Tree = Tree->brother; continue;}
            (*ngroot)++;
            uint64_t * ptu = (uint64_t*) malloc(nbcases*SUL);
            creabitfieldr(tindices,mod,Tree->link_rn->tab_data,ptu,nbcases-1,reste);
            Tree->link_rn->sup = sum;
            Tree->link_rn->tab = ptu;
            Tree = Tree->brother;
          }
        nb_freq += (*ngroot); 
      }
}

/// Regler Probleme avec l elagage 

void depthwalk (pnodes * Tree, uint64_t * cand)
{ 
  bool d = !Tree->brother;
  if (!Tree->link_rn->sup) {if (!d) {depthwalk(Tree->brother, cand);} return;}
  int supcand = 0;

  auto * tib = Tree->link_rn->tab;
  auto * ed = &cand[nbc];

  for (auto * tob = cand; tob != ed; ++tob,++tib )
    {supcand += __builtin_popcountl(*tob&*tib);}
 
  if (supcand<minsup) 
  {Tree->link_rn->sup = 0;
  if (!d) depthwalk(Tree->brother, cand);
  Tree->link_rn->sup = 1;
  return;
  }
           
  nb_freq++;
  if (!Tree->son)
    {Tree->son = (pnodes*)malloc(SPN);
     Tree->son->sup = supcand;
     Tree->son->link_rn = curlitem;
     Tree->son->son = NULL;
     Tree->son->brother = NULL;
     if (!d) {depthwalk(Tree->brother, cand);}
     return;
    }   
  pnodes * trpnodes = Tree->son;
  uint64_t  * ptab = (uint64_t*) malloc (nbc*sizeof(uint64_t)) ;
  auto * pptab = &ptab[0];
  tib = Tree->link_rn->tab;
  for (auto * tob = cand; tob != ed; ++tob,++tib,++pptab ){*pptab = *tob&*tib;}        
  depthwalk(trpnodes,ptab);
  free(ptab);
  Tree->son = (pnodes*)malloc(SPN);
  Tree->son->sup = supcand;
  Tree->son->link_rn = curlitem;
  Tree->son->brother = trpnodes;
  Tree->son->son = NULL; 
  if (!d) depthwalk(Tree->brother, cand);
  return;       
}


void rootwalk (pnodes * root )
{ pnodes * temp = root->brother;
  pnodes * trpn = NULL;

  while (temp)
    { 
      if (temp->link_rn->sup)
        { 
          trpn = temp->son;
          if (trpn == NULL)
            {
              temp->son = (pnodes*)malloc(SPN);
              temp->son->sup = temp->link_rn->sup;
              temp->son->link_rn = curlitem;
              temp->son->son = NULL;
              temp->son->brother = NULL;
              temp = temp->brother;
              continue;
            }
          depthwalk (trpn,temp->link_rn->tab);
          temp->son = (pnodes*)malloc(SPN);
          temp->son->sup = temp->link_rn->sup;
          temp->son->link_rn = curlitem;
          temp->son->son = NULL;
          temp->son->brother = trpn;
          temp = temp->brother;
          continue;
        }
      temp = temp->brother;

    }

}



double Bitprefrec (uint64_t** Bitdata, std::vector<std::string>& varnames, int maxul, int nvar, char ordre, pnodes ** ret_alpha, rnodes ** rnodes_ret ,int nb_var_moving, int nb_iter, char order,int & nadd, int & ndel)
{  curlitem = 0; nb_freq = 0;
 int sum = 0;
 auto a = 0; int b= 0;
 auto e = 0;
 auto n = 0;
 uint64_t ** Bdata = Bitdata ;
 auto * pul = &Bitdata[0][0];
 auto * pulb = &Bitdata[0][0];
 auto * pulc = &Bitdata[0][0];
 std::cout <<"Done The Supportvalue is " << minsup << ". Start PrefRec-Mooving with  " << std::endl;
 
 (*ret_alpha) = (pnodes*)malloc(sizeof(pnodes));
 pnodes * rootalpha = (*ret_alpha);
 rootalpha->brother = NULL;
 rootalpha->son = NULL;
 pnodes * lastree = rootalpha;
 pnodes * ttree = NULL;



  (*rnodes_ret) = (rnodes*)malloc(sizeof(rnodes));
 rnodes * rnodes_alpha = (*rnodes_ret); 
 rnodes_alpha->sup_data = 0;
 rnodes_alpha->brother = (rnodes*)malloc(sizeof(rnodes));
 rnodes_alpha->brother->sup_data = maxul*64 +100;
 rnodes_alpha->brother->brother = NULL;
 lastree->link_rn = rnodes_alpha;

 rnodes * trn = NULL;

 int j, ngroot;
 int m = 0;
 uint64_t u; 

 while (m < nb_var_moving)
  { 
    
    nb_freq++;   
    pul = Bdata[j];
    sum = 0;
    e =0;
    for (a = 0; a < maxul;a++)
      {sum+= __builtin_popcountl(pul[a]);}
    if (sum < minsup) {j++; continue;}
    m++;
 
    int *  indices = (int*)malloc(sum * SI);
    int * mod = (int*)malloc(sum * SI);
    for (a=0; a < maxul; a++)
      { u = pul[a];
        for (n = 0; n < 64; n++)
          { if ((u >> n) & 1UL) {indices[e]= a; mod[e++]=n;} }
      }

      
    rnodes * root_ritem [j+1];    
    nbc = sum/64;
    int reste = sum-64*nbc;
    if (reste) nbc ++;
    
    int s = 0;
    ttree = rootalpha->brother;
    while (ttree)
      {  
       s = 0; pulb = ttree->link_rn->tab_data;
       for (n = 0; n < maxul; n++) { s += __builtin_popcountl(pul[n]&pulb[n]);}
       ttree->link_rn->sup = s;
       ttree = ttree->brother;
      }  

    ngroot= 0;
    
    init_Rnodes(rootalpha,&ngroot,indices,mod,nbc,reste);
    free(indices);
    free(mod);
    
    rnodes * rnt = (rnodes*)malloc(sizeof(rnodes));
    rnt->sup_data = sum;
    rnt->tab = NULL;
    rnt->brother = NULL;
    rnt->tab_data = Bdata[j];
    rnt->var = &varnames[j];
    rnt->cor_bsb = NULL;

    

    curlitem = rnt;
    pnodes * nouvnod = (pnodes*)malloc(SPN);
    nouvnod->brother=rootalpha->brother;
    nouvnod->son=NULL;
    nouvnod->link_rn = rnt;
    nouvnod->sup = sum;

    j++;

    lastree->link_rn->cor_bsb = nouvnod ;

  
    set_rnodes_in_order (rnodes_alpha,rnt);
    rootwalk(rootalpha);
    ttree = rootalpha->brother;
    while (ttree)
      {
        if (ttree->link_rn->tab) {free(ttree->link_rn->tab); ttree->link_rn->tab = NULL;}
        ttree = ttree->brother;
      }
        
     rootalpha->brother = nouvnod;
     lastree = nouvnod;
    }



  int nb_del = 0;
  int maxiter = j+nb_iter;
  int n_freq_step = nb_freq;
  std::cout << "End of the first step, now starting the mooving application \n";
  nb_freq = 0;
  auto start = std::chrono::system_clock::now();
  for (j; j < maxiter;j++)
    { 

    nb_freq++;   
    pul = Bdata[j];
    sum = 0;
    e =0;
    for (a = 0; a < maxul;a++)
      {sum+= __builtin_popcountl(pul[a]);}
    if (sum < minsup){j++; continue;}
 
    int *  indices = (int*)malloc(sum * SI);
    int * mod = (int*)malloc(sum * SI);
    for (a=0; a < maxul; a++)
      { u = pul[a];
        for (n = 0; n < 64; n++)
          { if ((u >> n) & 1UL) {indices[e]= a; mod[e++]=n;} }
      }


    rnodes * root_ritem [j+1];    
    nbc = sum/64;
    int reste = sum-64*nbc;
    if (reste) nbc ++;
    
    int s = 0;
    ttree = rootalpha->brother;
    while (ttree)
      {  
       s = 0; pulb = ttree->link_rn->tab_data;
       for (n = 0; n < maxul; n++) { s += __builtin_popcountl(pul[n]&pulb[n]);}
       ttree->link_rn->sup = s;
       ttree = ttree->brother;
      }  

     
    ngroot= 0;
    init_Rnodes(rootalpha,&ngroot,indices,mod,nbc,reste);
    free(indices);
    free(mod);
    rnodes * rnt = (rnodes*)malloc(sizeof(rnodes));
    rnt->sup_data = sum;
    rnt->tab = NULL;
    rnt->brother = NULL;
    rnt->tab_data = Bdata[j];
    rnt->var = &varnames[j];
    rnt->cor_bsb = NULL;
    curlitem = rnt;
    
    pnodes * nouvnod = (pnodes*)malloc(SPN);
    nouvnod->brother=rootalpha->brother;
    nouvnod->son=NULL;
    nouvnod->link_rn = rnt;
    nouvnod->sup = sum;
    lastree->link_rn->cor_bsb = nouvnod ;
  
  set_rnodes_in_order (rnodes_alpha,rnt);
  
  rootwalk(rootalpha);
  

  ttree = rootalpha->brother;
  while (ttree)
    {
      if (ttree->link_rn->tab) {free(ttree->link_rn->tab); ttree->link_rn->tab = NULL;}
      ttree = ttree->brother;
    } 
    
    rnodes * todel = rnodes_alpha->brother;
    
    if (todel!= rnt)
      { 
        rootalpha->brother = nouvnod;
        lastree = nouvnod;
        
        erase_var (todel, &nb_del);
        
        rnodes_alpha->brother = todel->brother;
        free (todel);
         continue;
      }
    erase_nouvar (rootalpha,todel,&nb_del);  
    rnodes_alpha->brother = todel->brother;

    free (todel);
    free(nouvnod);
    nb_del++;
  }  

  ndel = nb_del;
  nadd = nb_freq;
  std::chrono::duration<double> diffe= std::chrono::system_clock::now() - start;
  std::cout << "End PrefRec-Mooving,   " << nb_freq <<" frequents have been added and <<" << nb_del << " deleted in " <<diffe.count() << " secs \n";
  std::cout << "it stays " << n_freq_step - nb_del + nb_freq << " frequent sets in the Tree \n ";

  return diffe.count();

}