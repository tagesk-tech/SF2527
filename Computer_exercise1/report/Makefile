JOBNAME := SF2527_CE1_Report

.PHONY: pdf clean

pdf:
	mkdir -p build
	latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build -jobname=$(JOBNAME) main.tex
	cp build/$(JOBNAME).pdf $(JOBNAME).pdf

clean:
	latexmk -C -outdir=build -jobname=$(JOBNAME) main.tex
