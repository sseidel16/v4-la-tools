# v4_analysis

This software can perform analysis or generate locating arrays.

All commands are from the bash shell within the main directory.

COMPILE:
A Makefile exists for compilation. Simply type:
$ make
If any errors occur, you can clean with:
$ make clean

ANALYSIS:
The analysis part of the software requires a valid locating array file, factor data file, and responses directory.
You can execute it by typing:
$ ./Search [LocatingArray.tsv] [FactorData.tsv] analysis [ResponsesDirectory] [response_column] [1/0 - perform log on responses] [nTerms] [nModels] [nNewModels]

FIXLA: (NOT tested extensively with groupings but SHOULD work)
This section fixes broken locating arrays by successively adding more rows.
A score of the array is kept (lower is better, 0 is valid locating array).
The software will print lower and lower scores until it reaches 0 and the valid locating array will be written to a file.
Pros: creates a LA by definition
Cons: the LAs may barely be valid
$ ./Search [LocatingArray.tsv] [FactorData.tsv] fixla [FixedOutputLA.tsv]

MTFIXLA:
This section fixes broken locating arrays by adding chunks of rows and randomizing them.
The randomization follows the Moser-Tardos idea.
A special variable k specifies how many differences are required for every pair of columns in the CS matrix.
This then produces stronger LAs.

REORDERROWSLA:
Reorder rows in LA by contributions.
Rows at the top of the LA are the ones contributing most.
Rows at the bottom often have 0 contribution and can be completely removed.
Under development...

Examples:
The following performs analysis on the LARGE simulated dataset (13 factors per model) (50 and 50 works pretty well):
./Search LA_LARGE.tsv Factors_LARGE.tsv analysis responses_LARGE Throughput 1 13 50 50
The following creates a locating array for the LARGE simulated dataset:
$ ./Search LA_LARGE_HEADER.tsv Factors_LARGE.tsv fixla OUT.tsv
