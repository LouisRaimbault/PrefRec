## Dataset test 

The purpose of this folder is to access various databases, whether synthetically generated or not.

If you want more information about the simulations results please have a look to the simulation section part of : (1)
The databases are in transaction format as defined below.

|transactions|
|------------|
|a,b|
|c,d|
|a,c|
|a,b,d|
|c|


They are currently in a compressed .data.gz format.
If you are using linux, you can use the command "$ gunzip 'file' " to get the classic data file.


### Kind of Database 


There are 3 main types of bases, synthetic independent and dependent as well as some known databases, reprocessed to be easily tested
The methods for generating synthetic bases are defined in (1). The table below gives the indications to better understand the names of the databases.
The ARn_4000_100_me658.txt and ARn_6000_100_me658.txt are separated into 2 files 1 and 2, each comprising half of the simulated transactions. Once unzipped you can get the full base with the command "cat file one file two > destination"



### Databases name indicators :
|indic|note|
|--------------------|--------|
|    B"format"   | format can be "d" for dependant, "i" for independant  |    
|    n"p" | p is the number of all different variables you can find in the transactions dataset|   
|    N"l"  | l in the number of observations/transactions divided by 1000 of the dataset | 
|    Me"m" | m is the average number of item per transaction|   
