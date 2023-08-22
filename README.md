## Haiti Article

This repository contains the necessary code to reproduce the article: "Informing policy via dynamic models: Cholera in Haiti". 

### TODO: 

- [ ] **DON'T FORGET**: Before re-running everything, remove `ncores` argument from `fit_haiti3` and `fit_haiti1`. This argument doesn't do anything, but was left in so as to not trigger the `bake` recalculation of the file, which would be inconvinient when working on the manuscript. 
- [ ] Make sure supplement is up to date 
   - [ ] Fix the supporting-information section in manuscript 
   - [ ] write an acknowledgement section 
   - [ ] Add model diagrams to supplement, and mention them in the manuscript. 
- [ ] Fix author affiliations and emails. 
- [ ] Double check all model equations match the code 
- [ ] Write more about why model 3 is a poor description of the data, even if this is just added to the supplement. 
- [ ] Update Model 3 initialization model description. 
- [ ] Final things: 
   - [ ] Check all journal requirements are met
   - [ ] Final proof-read.
   - [ ] Delete all commented out text and code 
   - [ ] Re-run everything and submit. 


- [x] Remove parallel back-ends from inside functions in `haitipkg`. 
- [x] Move all images and tables to the correct spot (directly after referencing them).
- [x] Have to reference all tables. When should we reference Table 2 (parameter table). How about the table of likelihoods? 
- [x] Remove vaccination scenario projection figures (maybe move into the supplement later).
- [x] Check we meet following requirement: "Please refrain from using math mode for non-math content, such as chemical formulas. For example, please use "CO\textsubscript{2}" instead of "$\mathrm{CO}_2$"
- [x] Fix the formatting of the parameter table. 
- [x] Write an author summary for the paper 
- [x] Description of theta (vaccine efficacy) in model 2 is missing. It's only in the equation, not even in the table. Make a mention of where this can be found.
- [x] Check that everything in the supplement material is referenced in the ms
- [x] Code has been refactored to calculate in the manuscript the likelihood of Lee et al model 1. Fix the supplemental material so that it doesn't have to redo these computations, then remove redundant comments from ms.Rnw
- [x] Check that all equations end in comma or period, remove unnecessary indentation after equations, and capitalize following text accordingly. 
- [x] Check that all in text citations are refer to the authors rather than the article (eg: in lee20 -> by lee20). 
- [x] Convert urls to xurls in .bib, except for packages. 
- [x] Capitalize title of all articles in .bib.
