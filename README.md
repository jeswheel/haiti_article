## Haiti Article

This repository contains the necessary code to reproduce the article: "Informing policy via dynamic models: Cholera in Haiti". 

The `ms.Rnw` file took 4 days, 9 hours, and 2 minutes to compile using 36 cores. It required 42 GB of memory. 

The `si/si.Rnw` file took 2 days, 11 hours, and 5 minutes to compile using 36 cores. It also required 36 GB of memory. 
In order for this file to compile, it is assumed that the data files created from `ms.Rnw` already exist. 

### How to recreate documents

First, the `ms.Rnw` needs to be knit into a `.tex` file. 
This is where the actually `R` code is run, and is very computationally intensive. 
As noted above, on 36 cores, this took longer than 4 days. 
Knitting the file will create all of the `R` data files that will then be loaded 
later. 
These files are used to create tables, figures, and output numbers in the 
manuscript. 
To avoid needing to recreate the `R` data files every time the code is ran, 
we use the function `pomp::bake`; this function will check if the code needs to 
be re-ran by checking if the code has changed or the reproducible seed before
the `pomp::bake` call has changed. 
If there are no changes, and the data file already exists, then the function 
simply loads existing results rather than recomputing them. 
Therefore if there are no changes to the `R` code that would force a recompute 
by `pomp::bake`, and the computed `R` data files already exist, then it only 
takes a few seconds to knit the document. 

### Plos Comp Bio submission

Once the document has been knit into a `.tex` file, several steps are taken to 
prepare the document for submission to Plos Comp Bio. 
For example, no graphics should be included in the submission, and all figures
need to be Tiffs. 
Also, the `.bbl` file should be included in the `.tex` document and not be 
loaded later. 

The make file contains a command `submission` which simplifies these steps. 
Running `make submission` will create a sub-folder called `submission/` that 
contains the `ms-submission.tex` file with all of the requirements. 
It also creates the corresponding `ms-submission.pdf` so that you can check to 
make sure that everything is working as expected. 
A final step that is not currently covered by `make submission` is to convert 
the parameter table into a figure.
This is necessary because the table contains some nested tabular environments 
(created by `\multirow`). 
To do this, I simply copy `submission/ms-submission.tex`, remove fancy header /
foot and page numbers, remove the caption from the table, compile into a pdf, 
and save the single page with the parameter table to it's own `.pdf` file 
named `submission/paramTab.pdf`. 
I can then convert this to a `.tiff` by running `make submission/paramTab.tiff`. 

Finally, it is necessary to manually change the table -> figure in the 
`submission/ms-submission/tex` file, and reference the table as a figure. 

### ArXiv submission

This is primarily documented for personal benefit. The `make ArXiv` command will 
create a `.zip` file that should be able to be submitted to ArXiv. 
For this to work, you need the `.bbl` files for both `ms.tex` and `si/si.tex`, 
so these should be compiled, but not cleaned, prior to running `make ArXiv`. 
A more advanced modification to the make file would make sure these exist, but 
because this is just for personal use, it felt unecessary to spend the additional 
effort to create the make file in this case. 

### TODO: 

Nothing to do! 
