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

**PrefRules:**
* (recursiv) Mining confident rules with a minConf value
* 7 differents possibles kind of supplementary coefficient
* Possibility of adding additional constraints to the extraction of the rules, based on the proposed coefficients
* Possibility to select target items or itemset


**PrefRules coefficient :** 
* Conf2 (usually know as minconf)
* Conf1
* Power2
* Power1
* Covariance
* Correlation
* Kappa
* Maxwell-Pilliner 

**Prefrules global coefficient :**
* Mean = (1/4) (conf1+conf2+power1+power2)
* Strong = 1- sup(antecedant)-sup(consequent) + 2* sup(set)

## Creating binaries and getting started
```
cd Prefrec_src && make
cd Prefrules_src && make
./Prefrec <transaction dataset> < d=item_delimitator> <s=minimal_relative_support "s"> <output_file_set_infos> <output file for coefficient>
./Prefrules <file set coefficient> <file set coefficient rules> <number of n condition coefficient> <coeff1=value> ... <target item> <path for output> <coeffn=value> ... <number of m extract coefficient> <coeff1> ... <coeffm> <bool ok size> <bool ok sup> <bool ok global_indic>

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

The section https://github.com/LouisRaimbault/PrefRules_Strat/tree/main/databases_test presents some larger databases used for efficiency simulations. If you want more information, you can read the associated read me.


### Example
```
./Prefrec ../sample/input/Fruits.txt d=, s=0.20  \\ Do the extraction on the bases "fruits" with a relativ minSup of 0.20, item are delimited by a ",". 
./Prefrec ../sample/input/Fruits.txt d=, s=0.20 ../sample/output_Prefrec/infoset \\ Do the extraction and write frequent set informations in infoset   
./Prefrec ../sample/input/Fruits.txt d=, s=0.20 ../sample/output_Prefrec/infoset ../sample/output_Prefrec/genrules \\ Do the extraction, write frequent set informations in info set,  and create 2 files , genrules.txt and genrules_item.txt, necessary to use Prefrules

./Prefrules ../sample/output_Prefrec/genrules.txt ../sample/output_Prefrec/genrules_variables.txt 2 Conf2=0.5 Conf1=0.2 all_tgts \\ Calculate all confident rules with a minconf of 0.5 , then valid rules with conf1>=0.2
./Prefrules ../sample/output_Prefrec/genrules.txt ../sample/output_Prefrec/genrules_variables.txt 2 Conf2=0.5 Conf1=0.2 apple, \\ Calculate confident rules with a minconf of 0.5 and having an antecedant or consequebt containing apple, then valid rules with conf1>=0.2
./Prefrules ../sample/output_Prefrec/genrules.txt ../sample/output_Prefrec/genrules_variables.txt 1 Conf2=0.5 all_tgts ../sample/output_Prefrules/inforules 2 Conf2 Power2 1 1 0 \\ Calculate all confident rules with a minconf of 0.5 and write confident rules with their Conf2,Power2 and their size and supports values in inforules
./Prefrules ../sample/output_Prefrec/genrules.txt ../sample/output_Prefrec/genrules_variables.txt 2 Conf2=0.5 Conf1=0.2 1 all_indicators apple, 1 1 1 ../sample/output_Prefrules/inforules  \\  Calculate confident rules with set containing the item "apple" with a minconf of 0.75 and write all the rules and the 8 indicator values with their size, their support and the 2 global indicators
``` 


## Parameters for Prefrec :
|param|required|note|
|--------------------|--------|--------|
|    transaction dataset    |    yes    | The Dataset transaction  |  
|    item delimitator "d="   |    yes    | the item separator you use with your dataset transaction | 
|    minimal relative support "s="   |    yes    | the minimal relativ support you wish     | 
|    ordering frequent 1-itemset "o="   |    no    | n for unordered, i for ascendant order, d for decreasing order      | 
|    output file set infos    |    no    |  For file frequent set informations, put the directory file (with no extention type )    | 
|    output file for coefficient    |    no    |  If you want to use Prefrules, put the directory file (with no exention type)| 



## Parameters for Prefrules :
|param|required|note|
|--------------------|--------|--------|
|    file set coefficient   |    yes    | The genrules_file.txt created with Prefrec   |
|    file set coefficient item |    yes    |  The genrules_file_item.txt created with Prefrec  | 
|    number of coefficient condition |    yes    |  The number of coefficient  |
|    coefficient condition |    yes    |  all the coefficient with their value , ex Conf1=0.1  |
|    (1)output file coefficient infos    |    no    |  If you want the output file rules informations, put the directory file (with no exention type)    |
|    number of extracted coefficient |    yes if (1)   |  The number of extracted coefficient, if you wan't all type 1  | 
|    extracted coefficient |    yes if (1)    |  all the extracted coefficient, if you wan't all type "all_indicators"  |
|    target +','  |    yes if (1)   | the target item. if you wan't to set a particular target, otherwise type "all_tgts"  |
|    ok sup |    yes if (1)   |  if you wan't the support values for antecedant consequent and set  |
|    ok size |    yes if (1)   |  if you wan't the size values for antecedant consequent and set  |
|    ok global indic |    yes if (1)   |  if you wan't the 2 globals indic, Mean and Strong  |
  





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




