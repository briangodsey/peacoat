About PEACOAT
-------------

[This paper in the Journal of Integrative
Bioinformatics](http://journal.imbio.de/article.php?aid=227) presents
a Prediction-Enhanced Algorithm for Clustered, Ordered Assessment of
Targeting in miR-mRNA interactions called *PEACOAT*. Please read the
paper for details about the algorithm.


Using PEACOAT
-------------

> **WARNING**: This code is messy and WILL be difficult to use. I am a
> mathematician, not a software developer, and I have put significant
> time into making this code work as well as possible. However, I am
> sure that it won't work the first time for everyone who tries it. If
> you are having difficulties after reading this documentation as well
> as the referenced description files, feel free to contact me with
> details about the problems you are having. Likewise, if you would like
> to contribute to *PEACOAT* and help make it better, I would greatly
> appreciate it.

This *Git* repository contains code for *PEACOAT* that can reproduce
the results in the above paper on *R* version 3.0. Many packages are
required, but for now I won't list them here; they will be obvious
upon running the code.

To reproduce the results on the Myeloma data referenced in the paper,
use the R script

````
./myeloma/script-myeloma.R
````

To simulate results as in the paper, use the R script

````
./simulation/script-simulations.R
````

with the parameter **geo.files** set to the appropriate paths to the
dowloaded files from
[GSE17498](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17498).

Note, however, that one should be careful using these scripts, because

* they can take quite a long time to run
* there are many data file dependencies

So, please be sure to read and understand the code before running the
scripts. I don't recommend running the entire script at once until you
have run the code in sections at least once.


User-defined parameters
-----------------------

The file

````
./myeloma/README_userDefinedParameters_myeloma.R
````

and its counterpart in the simulation/ folder contain all of the
user-defined parameters contains the configuration settings that a
typical user might need. This is a fairly complex file and attention
is needed to select the proper parameter values. A lot of comments
have been included in the file for this reason.

Parameters must be set according to the expression data type that will
be used. The following data types have code that supports them:

* GEO
* expression matrix
* Agilent
* Affymetrix
* simulation

Each may require certain *R* packages.

Among the many user-defined parameters, there are two that warrant
special mention:

* typefollows.generator
* typeadder.script

Each of these is a path to a filename that handles the *types* within
the expression data. For example, a data set might have three stages:

* embryo
* youth
* adult

and we would like to consider the stages in their natural order. To
accomplish this, we would either create the variable *typefollows* to
reflect this order (see the parameter file for examples) or we would
write a script that generates the variable. The parameter
*typefollows.generator* should be set to the name of this
script. Often in simple cases, a script is not necessary, but in more
complex cases it can be helpful. Likewise, the parameter
*typeadder.script* should be set to the name of a script that adds the
types to the expression object (in the *targets* slot), and is
necessary only if the types are not present already.

See [my previous paper in PLoS
ONE](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0051480)
for more info about types and orders.


File dependencies
-----------------

Currently, starting from the repo directory---i.e. where the .git/
folder is located---I am using the following directories outside of
the repo:

````
../databasefiles
../geo-data
../generatedfiles
````

The **databasefiles/** contains

````
predicted/
    miRanda/
        hg19_predictions_S_C_aug2010.txt
        mm9_predictions_S_C_aug2010.txt
    targetscan/
        miR_Family_Info.txt
        human/
            Conserved_Site_Context_Scores.txt
        mouse/
            Conserved_Site_Context_Scores.txt

validated/
     miRecords_version3.csv
     miRWalk_validatedTargets.txt
     TarBase_V5.0.csv
````

where the files of predicted miR targets can be downloaded from the
[miRanda downloads
page](http://www.microrna.org/microrna/getDownloads.do) and the
[TargetScan](http://www.targetscan.org) downloads page for the
organism in question.

The locations of the *predicted* files can be modified in the file

````
./src/loadDatabasePredictions.R
````

while the paths to the *validated* files as well as
*miR_Family_Info.txt* can be nodified in

````
./src/modelDataInitializations.R
````

Generated files
---------------

The folder ../generatedfiles/, whose location can be changed in
./src/modelDataInitializations.R, is where intermediate,
code-generated files are to be stored. Mainly, this is because
processing the database prediction files and other initializations
take some time, so the intermediate results are stored to save time
on future script runs.

> Note that the script checks for the existence of certain intermediate
> files, uses them if they are present, and generates them if they are
> not. Thus, to re-calculate a certain intermediate file, either delete
> it or re-name it.


Questions, concerns, and collaborators
--------------------------------------

Certainly I've skipped many details about the code here, but please
feel free to contact me with questions or suggestions. I would also
appreciate contributions, pull requests, etc.
