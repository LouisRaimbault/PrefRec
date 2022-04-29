#include "Prefrec_init.h"
#include "Prefrec_functions.h"

#pragma optimization_level 3
#pragma GCC optimize("Ofast,no-stack-protector,unroll-loops,fast-math,O3")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx")
#pragma GCC optimize("Ofast")
#pragma GCC target("avx,avx2,fma")
#pragma GCC optimization ("unroll-loops")

int main (int argc , char ** argv )
{
  nb_freq = 0;
  char  * transacpath = (char*)argv[1];
  std::string delimparam ((char*)argv[2]);
  std::string relsupparam ((char*)argv[3]);
  int nb_var_moving = atoi((char*)argv[4]);
  int nb_iter = atoi((char*)argv[5]);
  std::string orderparam ;
  
  char delim ;
  char ordre = 'u';
  for (int i = 0; i < 64; i++) {PT_ULTAB[i] = (1UL <<i);}
  if (delimparam[0]=='d' && delimparam[1]=='=')
    {
      delim = delimparam[2];
    }
  else { std::cout << "Wrong second argument, please set d=yourdelim for separate item transactions" << std::endl; return 0;}

  if (relsupparam[0]=='s' && relsupparam[1]=='=')
    {relsupparam.erase(relsupparam.begin(),relsupparam.begin()+2);}
  else {std::cout << "Wrong third argument, please set s=yourrelativeminSup" << std::endl; return 0;}
  short n_param = 6;      
  if (argc > 6 ) 
    {orderparam = std::string ((char*)argv[6]);
     if (orderparam[0] == 'o' && orderparam[1] =='=') {ordre = orderparam[2]; n_param++; }
    } 
  int nrows=0; int nvar = 0; int maxul = 0;
  uint64_t ** Bitdata = NULL;
  std::vector<std::string> varnames;
  Init_data(transacpath,nrows,maxul,nvar,Bitdata,varnames,delim);
   double relativeminsup = std::stod(relsupparam);
   minsup = nrows* relativeminsup;
   pnodes * tab_root_pnodes  = NULL;
   rnodes * rnodes_alpha = NULL;

   int ndel = 0;
   int nadd = 0;
   double time_n = Bitprefrec(Bitdata, varnames,maxul,nvar,ordre,&tab_root_pnodes,&rnodes_alpha,nb_var_moving,nb_iter, ordre, nadd,ndel);

   for (int s =0; s < nvar;s++)
   { free(Bitdata[s]);}
     free(Bitdata);
    
   
   if (nb_freq == 0) { std::cout << "Not a single frequent itemSet, please try with a lower relative minSup value. " << std::endl; return 0;}
  if (argc <= n_param) 
    { std::cout << "Erase Frequent set .. \n " ; 
      erase_freq (tab_root_pnodes); 
      rnodes * tr_rn = rnodes_alpha;
      while (tr_rn)
        {
          tr_rn = rnodes_alpha->brother;
          free(rnodes_alpha);
          rnodes_alpha = tr_rn;
        }
      std::cout << "done" << std::endl;
      return 0;
    }  
  
  if (argc == n_param+1)
    { 
      double nrs = (double)nrows;
      char  * argoutpath = (char*)argv[n_param];
      std::string outpath (argoutpath);
      outpath += ".txt";
      std::cout << "Filling Frequent sets informations in " << outpath << " and delete Tree ... \n";
      std::ofstream flux_set;
      flux_set.open(std::string(outpath));
      flux_set << "nameset\tSupport\trelative_support\tsize\n";
      int size = 1;
      std::string ststring = "";
      extract_and_erase_freq (tab_root_pnodes->brother,ststring,&varnames[0],nrs,size,flux_set);
       rnodes * tr_rn = rnodes_alpha;
      while (tr_rn)
        {
          tr_rn = rnodes_alpha->brother;
          free(rnodes_alpha);
          rnodes_alpha = tr_rn;
        }
      std::cout << "done \n";
      flux_set.close();  tr_rn = rnodes_alpha;
  while (tr_rn)
    {
      tr_rn = rnodes_alpha->brother;
      free(rnodes_alpha);
      rnodes_alpha = tr_rn;
    }
      return 0;
    }

  int m = 0;
  int seuil_max = maxul*64+80;
  rnodes * tr_rn = rnodes_alpha->brother; 
  while (tr_rn->sup < seuil_max) {varnames[m] = *tr_rn->var; tr_rn->lab = m; m++;}
  double nrs = (double)nrows;  
  char  * argoutpath = (char*)argv[n_param];
  n_param++;
  char  * argoutpathrules = (char*)argv[n_param];
  std::string outpath (argoutpath);
  std::string outpathrules (argoutpathrules);
  outpath += ".txt";
  outpathrules += "_variables.txt";
  std::ofstream flux_set;
  std::ofstream flux_rules;
  flux_rules.open(outpathrules);
  std::cout << "Writing frequent 1-item name in " << outpathrules << "\n";
  for(int j = 0; j < m;j++){ flux_rules << varnames[j] <<"\n";}
  flux_rules.close();

  outpathrules = std::string(argoutpathrules);
  outpathrules += ".txt";
  std::cout << "Filling Frequent sets informations in " << outpath <<  " informations for PrefRules in " << outpathrules <<" and delete Tree ... \n";
  flux_set.open(outpath);
  flux_rules.open(outpathrules);
  flux_set << "nameset\tSupport\trelative_support\tsize\n";
  int size = 1;
  std::string ststring = "";
  extract_and_erase_freq (tab_root_pnodes->brother,ststring,&varnames[0],nrs,size,flux_set,flux_rules);
  std::cout << "done \n";
  tr_rn = rnodes_alpha;
  while (tr_rn)
    {
      tr_rn = rnodes_alpha->brother;
      free(rnodes_alpha);
      rnodes_alpha = tr_rn;
    }
  flux_set.close();
  flux_rules.close();

  return 0;
}
