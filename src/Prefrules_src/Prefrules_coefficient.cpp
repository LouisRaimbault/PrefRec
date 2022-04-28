
#include "Prefrules_coefficient.h"

double conftwo (double  freq_ante, double  freq_complem, double  freq_set )
{
  return freq_set/freq_ante;
  
}

double powertwo (double  freq_ante, double  freq_complem, double  freq_set )
{
  if (abs(freq_ante -1) < 0.0000001) {return 1;}
  return (1-freq_ante - freq_complem + freq_set)/(1-freq_ante);
}

double confone (double  freq_ante, double  freq_complem, double  freq_set )
{
  if (abs(freq_complem -1) < 0.0000001) {return 1;}
  return (1-freq_ante - freq_complem + freq_set)/(1-freq_complem);
}

double powerone (double  freq_ante, double  freq_complem, double  freq_set )
{ 
  return freq_set/freq_complem;
}

double cov (double  freq_ante, double  freq_complem, double  freq_set )
{  
  return freq_set - freq_ante*freq_complem;
}

double corr (double  freq_ante, double  freq_complem, double  freq_set)
{ 
  double cov = freq_set - freq_ante*freq_complem;
  return cov/std::pow(freq_ante*(1-freq_ante)*freq_complem*(1-freq_complem),0.5);
}


double kappa (double  freq_ante, double  freq_complem, double  freq_set)
{ 
  double cov = freq_set - freq_ante*freq_complem;
  return cov/(freq_ante*(1-freq_complem)+freq_complem*(1-freq_ante));
}

double maxwellP (double  freq_ante, double  freq_complem, double  freq_set)
{
  double cov = freq_set - freq_ante*freq_complem;
  return (2*cov)/(freq_ante*(1-freq_ante)+freq_complem*(1-freq_complem));
}

