# v4_analysis

This software can perform analysis or generate locating arrays.

All commands are from the bash shell within the main directory.

Locating array examples are included in the LA/ directory.
Locating array format used by the software has changed over time, and a version on the first line of the LA is required for the software to accept it.
The LA version used by the software can be found in the LocatingArray.h header file.
Not all LAs in the LA/ directory use the latest version format (most use earlier versions at the time of writing).
This software only ckecks the first line of a LA to protect from version mistakes, and then assumes any LA passed as an argument is a correctly formatted locating array.
If a wrongly formatted locating array is passed to the software, undefined behavior will occur; often the software simply crashes.
The software is NOT tolerant to any mistakes in LA format.

## COMPILE
A Makefile exists for compilation. Simply type:
```
$ make
```
If any errors occur, you can clean with:
```
$ make clean
```

## SEARCH
Usage: `./Search [LocatingArray.tsv] ([FactorData.tsv]) ...`
Compilation generates an executable that can be run using `./Search`.
The executable requires, in general, 2 files: a locating array (LA) file and a factor data (FD) file.
The LA file uses factor and level indices that can be decoded by the FD file.
In addition, the FD file indicates if factors have numeric levels (useful for constraints) and what the float values are for each level.
Optionally, an empty string ("") can be provided instead of a FD file, and the software will assign all factors and (non-numeric) levels generic names (F0, F1, L0, L1...).
The seed used for all random samples throughout execution is obtained from the current time and printed immediately after execution begins.

## ANALYSIS
Usage: `... analysis [ResponsesDirectory] [response_column] [1/0 - perform log on responses] [nTerms] [nModels] [nNewModels]`
The analysis part of the software requires a valid locating array file, factor data file, and responses directory.
Analysis runs at a high level in `Search.cpp`, but it calls functions in `Model.cpp`.
It uses least squares to fit models, specifically QR decomposition, and this is all coded into `Model.cpp`.
You can execute it by typing:
```
$ ./Search [LocatingArray.tsv] ([FactorData.tsv]) analysis [ResponsesDirectory] [response_column] [1/0 - perform log on responses] [nTerms] [nModels] [nNewModels]
```

## CHECKLA
Usage: `... checkla [k Separation] [c Minimum Count]`
This section checks a LA file for a required separation and a minimum count for every interaction.
It performs three main checks.
First is a path check that uses the usual scoring for k separation and c minimum count.
Second is a weird scoring method that was originally used with FIXLA.
Third is a brute force check that produces the exact same score as the path check (at least it should, and it has for every LA I have checked).
I have not removed brute force because I am not 100% sure path check does not have some edge case bug and I want to spot it easily if it pops up.
However, it is placed last because it takes a long time and execution can simply be terminated anytime that somebody does not want to wait that long.

## FIXLA
Usage: `... fixla [FixedOutputLA.tsv]`
This section fixes broken locating arrays by successively adding more rows using the inital greedy construction approach in my thesis.
A score of the array is kept (lower is better, 0 is valid locating array).
The software will print lower and lower scores until it reaches 0 and the valid locating array will be written to a file.
Pros: creates a LA by definition, the LA is smaller than from other approaches, it is constructed fairly quickly.
Cons: the LA may barely be valid, does not incorporate separation or a minimum count for each interaction, cannot be used with constraints.
```
$ ./Search [LocatingArray.tsv] ([FactorData.tsv]) fixla [FixedOutputLA.tsv]
```

## MTFIXLA
This section fixes broken locating arrays by adding chunks of rows and randomizing them.
The randomization follows the Moser-Tardos idea.
A special variable k specifies how many differences are required for every pair of columns in the CS matrix.
This then produces stronger LAs.

## REORDERROWSLA
Usage: `... reorderrowsla [k Separation] [c Minimum Count] [ReorderedOutputLA.tsv]`
The argument k indicates the separation for the LA, and c indicates the minimum number of times each interaction must be covered.
Reorder rows in LA by contributions.
A contribution in a row is either interactions that differ when they have not already differed k times, or an interaction that is covered whenit has not already been covered c times.
Rows at the top of the LA are the ones contributing most.
Rows at the bottom often have 0 contribution and can be completely removed.
The generated output LA does not remove any rows but simply reorders them, and the console indicates which rows have contributions.

## CNERT LA Commands
The first analysis command below is for MOS and the second is for Exposure.
A log is used with Exposure, but no log command is given to the analysis because the log has already been taken in the response file.
This differs from throughput in Abraham's LA where no log is taken in the response file so a log must be taken by the analysis software.
```
./Search LA/LA_SMALL.tsv FD/Factors_SMALL.tsv analysis RE/responses_SMALL MOS 0 11 50 50
./Search LA/LA_SMALL.tsv FD/Factors_SMALL.tsv analysis RE/responses_SMALL Exposure 0 11 50 50
```

## Examples
The following performs analysis on the LARGE simulated dataset (13 factors per model) (50 and 50 works pretty well).
```
./Search LA/LA_LARGE.tsv FD/Factors_LARGE.tsv analysis RE/responses_LARGE Throughput 1 13 50 50
```
The following creates a locating array for the LARGE simulated dataset.
```
$ ./Search LA/LA_LARGE_HEADER.tsv FD/Factors_LARGE.tsv fixla OUT.tsv
```
The following check a LA that contains constraint groups, and check a LA for the LARGE simulated dataset.
The first requires its corresponding FD file because it uses constraints and numeric levels for its factors.
The second does not required its corresponding FD file because no factors have numeric levels for constraints and the LA can be checked without knowing the names of factors and levels.
```
./Search LA/LA_LARGE.tsv FD/Factors_LARGE.tsv analysis RE/responses_LARGE Throughput 1 13 50 50
./Search LA/LA_LARGE.tsv "" checkla 1 3
```
