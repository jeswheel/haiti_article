SUPFILES:= {mod1diagram,mod2diagram,mod3diagram,modelDetails,paramUncertainty,translateTab,measurementModels}
COMPILEFILES:= calibrateMod3Out.tex ReplicateLee20Out.tex confidenceIntervalsOut.tex initialValuesOut.tex
# INPUTFILES:= inputs/header.tex inputs/mod1diagram.tex inputs/mod2diagram.tex inputs/mod3diagram.tex inputs/modelDetails.tex inputs/measurementModels.tex inputs/initialValues.tex inputs/paramUncertainty.tex

create:
	$(MAKE) si.Rnw
	$(MAKE) sifiles
	$(MAKE) clean

# create: ReplicateLee20.Rnw calibrateMod3.Rnw
# 	if ! [ -d "si_files" ]; then mkdir si_files; fi
# 	$(MAKE) si.Rnw
# 	$(MAKE) clean

fresh: ReplicateLee20.Rnw calibrateMod3.Rnw $(INPUTFILES)
	rm -f si.tex calibrateMod3*.tex ReplicateLee20*.tex confidenceIntervals*.tex initialValues*.tex
	rm -rf si_files

si.Rnw: | $(COMPILEFILES) inputs/translateTab.tex
	Rscript -e "library(knitr); knit2pdf('si.Rnw')"

%Out.tex: %.Rnw
	if ! [ -d "si_files" ]; then mkdir si_files; fi
	Rscript -e "library(knitr); knit2pdf('$*.Rnw')"
	mv $*.pdf si_files/$*.pdf
	Rscript --vanilla clean_tex.R $*.tex
	rm $*.tex

clean:
	rm -f *.bbl *.log *-concordance.tex *.gz *.toc *.aux *.blg

si_files/%.pdf: inputs/%.tex utilityFuns.R make_si_file.R
	if ! [ -d "si_files" ]; then mkdir si_files; fi
	Rscript --vanilla make_si_file.R inputs/$*.tex
	pdflatex $*.tex
	if [ -a $*.tmp ] ; \
	then \
		bibtex $*.aux ; \
		pdflatex $*.tex ; \
		pdflatex $*.tex ; \
		rm $*.tmp ; \
	fi;
	mv $*.pdf si_files/$*.pdf
	$(MAKE) clean
	rm $*.tex

sifiles:
	for f in $(SUPFILES); do (make si_files/$$f.pdf); done


