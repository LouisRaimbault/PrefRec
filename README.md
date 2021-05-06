# PrefRec

## Description 

The purpose of this Git page is to access the uses of the **PrefRec** and **PrefRules** C++ programs, which can only run by using the GCC compiler until the version 1.2.
These 2 algorithms were built around the notion of *Mining Frequent ItemSets and associations rules*. This data analysis method was first introduced by *Agrawal et al. 1993* for mining transaction databases.

**Prefrec** and **Prefrules** use functions that allow recursive application by variable as well as an update method (new items can be added as long as the individuals remain the same). However, the current version (1.1) performs the usual applications of this type of Data Mining. Because it is complicated to predict the various wishes of users.
Please, Let us know if you want to get a specialized version for your use case, and quote my work if you want to reuse this code for your personal needs.

More versions will arrive soon for more precise and diverse applications.
If you would like to better understand their use and methods of operation, please have a look to : (link comming soon) (1)

**Prefrec:**
* (recursiv) Mining frequent itemSets with a relative Support
* Get an output Tsv file for information about frequent itemSet
* Get 2 files necessary to launch Prefrules 

**Prefrules:**
* (recursiv) Mining confident rules with a minConf value
* 4 differents possibles kind of supplementary coefficient
* Get an output txt files for informations about confident rules and the coefficients you choose 

**Prefrules supplementary coefficient :** 
* Covariance
* Correlation
* Kappa
* Maxwell-Pilliner 

## Creating binaries and getting started
```
cd src && make
./Prefrec <transaction dataset> < d=item_delimitator> <s=minimal_relative_support "s"> <output_file_set_infos> <output file for coefficient>
./Prefrules <file set coefficient> <file set coefficient rules> <c=minconf> <n=nbcoeff> <coeff 1 > ... <coeff n>  <output file coefficient info>

```

## Dataset format 

For this version, the Dataset format accepted is a file of type transaction, such that the items are separated by a choosen separator
and the transaction by a new line.Please, now that the sep must be one of the ASCII Chart. 
The following tab is an example with sep=, item a b c and d, and five transactions:



|transactions|
|------------|
|a,b|
|c,d|
|a,c|
|a,b,d|
|c|


The folder sample contain a simple small dataset test. You can use it for the example below to see how the software works.

The section https://github.com/LouisRaimbault/PrefRec/tree/main/databases_test presents some larger databases used for efficiency simulations. If you want more information, you can read the associated read me.


### Example
```
./Prefrec ../sample/input/Fruits.txt d=, s=0.20  \\ Do the extraction on the bases "fruits" with a relativ minSup of 0.20, item are delimited by a ",". 
./Prefrec ../sample/input/Fruits.txt d=, s=0.20 ../sample/output_Prefrec/infoset \\ Do the extraction and write frequent set informations in infoset   
./Prefrec ../sample/input/Fruits.txt d=, s=0.20 ../sample/output_Prefrec/infoset ../sample/output_Prefrec/genrules \\ Do the extraction, write frequent set informations in info set,  and create 2 files , genrules.txt and genrules_item.txt, necessary to use Prefrules

./Prefrules ../sample/output_Prefrec/genrules.txt ../sample/output_Prefrec/genrules_variables.txt c=0.75 n=0 \\ Calculate confident rules with a minconf of 0.75
./Prefrules ../sample/output_Prefrec/genrules.txt ../sample/output_Prefrec/genrules_variables.txt c=0.75 n=0 ../sample/output_Prefrules/inforules \\ Calculate confident rules with a minconf of 0.75 and write confident rules info in inforules
./Prefrules ../sample/output_Prefrec/genrules.txt ../sample/output_Prefrec/genrules_variables.txt c=0.75 n=2 cov corr ../sample/output_Prefrules/inforules \\  Calculate confident rules with a minconf of 0.75, compute their coeffiient cov and corr and write  informations in inforules.txt
```


## Parameters for Prefrec :
|param|required|note|
|--------------------|--------|--------|
|    transaction dataset    |    yes    | The Dataset transaction  |  
|    item delimitator "d="   |    yes    | the item separator you use with your dataset transaction | 
|    minimal relative support "s="   |    yes    | the minimal relativ support you wish     | 
|    output file set infos    |    no    |  For file frequent set informations, put the directory file (with no extention type )    | 
|    output file for coefficient    |    no    |  If you want to use Prefrules, put the directory file (with no exention type)| 



## Parameters for Prefrules :
|param|required|note|
|--------------------|--------|--------|
|    file set coefficient   |    yes    | The genrules_file.txt created with Prefrec   |
|    file set coefficient item |    yes    |  The genrules_file_item.txt created with Prefrec  | 
|    minconf  "c="  |    yes    | the minconf you want to use for prunning, if you dont want to, put -1  |
|    nbcoeff "n=" |    yes    | the number of supplementary coefficient you want to compute |   
|    type of coefficient    |    no if n=0    | the supplementary coefficient you want to compute | 
|    output file coefficient infos    |    no    |  If you want the output file rules informations, put the directory file (with no exention type)    |  





## Definitions of frequent itemSets :

Let us remind you the 2 mains definitions of this data analys method

Let I = {a1, . . . , an} be a finite set of items. A transaction database is a set of transactions T =
{t1, . . . , tN } where each transaction ti ⊂ I, 1 ≤ i ≤ N, represents a nonempty
subset of items. An itemset A is a subset of I; A is a k-itemset if it contains
k items. The support of an itemset A is denoted as supp(A) and is defined
as the number of transactions which contain A. The relative support of A is
freq(A) = supp(A)/N. A is frequent if freq(A) ≥ σ where σ is a user-specified minimum relative support threshold, called minSup.


**Definitions of confident associations rules  :**
An association rule is an implication A ⇒ B where A and B are two itemsets. The support of a rule A ⇒ B is defined as sup(A ⇒ B) = sup(A∪B).
The confidence of a rule A ⇒ B is defined as conf(A ⇒ B) = supp(A ⇒B)/supp(A) = freq(A∪B)/freq(A).
Considering a treshold minconf, a rule such that A ⇒ B is confident if conf (A ⇒ B) > minconf.

**Others criterion**

We propose some other criterion:
- Covariance 
- Correlation
- Kappa coefficient 
- Maxwell-Pilliner coefficient




