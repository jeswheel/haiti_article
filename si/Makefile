SUPFILES:= {mod1diagram,mod2diagram,mod3diagram,modelDetails,paramUncertainty,translateTab,initialValues,measurementModels}
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
	rm -f si.tex calibrateMod3*.tex ReplicateLee20*.tex
	rm -rf si_files

si.Rnw: | calibrateMod3Out.tex ReplicateLee20Out.tex inputs/translateTab.tex
	Rscript -e "library(knitr); knit2pdf('si.Rnw')"

calibrateMod3Out.tex: calibrateMod3.Rnw
	if ! [ -d "si_files" ]; then mkdir si_files; fi
	Rscript -e "library(knitr); knit2pdf('calibrateMod3.Rnw')"
	mv calibrateMod3.pdf si_files/calibrateMod3.pdf
	Rscript --vanilla clean_tex.R calibrateMod3.tex
	rm calibrateMod3.tex

ReplicateLee20Out.tex: ReplicateLee20.Rnw
	if ! [ -d "si_files" ]; then mkdir si_files; fi
	Rscript -e "library(knitr); knit2pdf('ReplicateLee20.Rnw')"
	mv ReplicateLee20.pdf si_files/ReplicateLee20.pdf
	Rscript --vanilla clean_tex.R ReplicateLee20.tex
	rm ReplicateLee20.tex

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


