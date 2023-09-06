create: ReplicateLee20.Rnw calibrateMod3.Rnw inputs/header.tex inputs/mod1diagram.tex inputs/mod2diagram.tex inputs/mod3diagram.tex inputs/modelDetails.tex inputs/measurementModels.tex inputs/initialValues.tex calibrateMod3Out.tex ReplicateLee20Out.tex inputs/paramUncertainty.tex inputs/translateTab.tex
	if ! [ -d "si_files" ]; then mkdir si_files; fi
	$(MAKE) si.Rnw
	$(MAKE) clean

fresh: ReplicateLee20.Rnw calibrateMod3.Rnw inputs/header.tex inputs/mod1diagram.tex inputs/mod2diagram.tex inputs/mod3diagram.tex inputs/modelDetails.tex inputs/measurementModels.tex inputs/initialValues.tex inputs/paramUncertainty.tex inputs/translateTab.tex
	rm si.tex calibrateMod3*.tex ReplicateLee20*.tex
	$(MAKE) si.Rnw
	$(MAKE) clean

si.Rnw: ReplicateLee20.Rnw calibrateMod3.Rnw inputs/header.tex inputs/mod1diagram.tex inputs/mod2diagram.tex inputs/mod3diagram.tex inputs/modelDetails.tex inputs/measurementModels.tex inputs/initialValues.tex inputs/paramUncertainty.tex inputs/translateTab.tex | calibrateMod3Out.tex ReplicateLee20Out.tex
	Rscript -e "library(knitr); knit2pdf('si.Rnw')"

calibrateMod3Out.tex: calibrateMod3.Rnw
	Rscript -e "library(knitr); knit2pdf('calibrateMod3.Rnw')"
	mv calibrateMod3.pdf si_files/calibrateMod3.pdf
	Rscript --vanilla clean_tex.R calibrateMod3.tex
	rm calibrateMod3.tex

ReplicateLee20Out.tex: ReplicateLee20.Rnw
	Rscript -e "library(knitr); knit2pdf('ReplicateLee20.Rnw')"
	mv ReplicateLee20.pdf si_files/ReplicateLee20.pdf
	Rscript --vanilla clean_tex.R ReplicateLee20.tex
	rm ReplicateLee20.tex

clean:
	rm -f *.bbl *.log *-concordance.tex *.gz *.toc
