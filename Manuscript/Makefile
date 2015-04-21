all: ohyou.pdf

ohyou.pdf: ohyou.tex
	pdflatex ohyou
	bibtex ohyou
	pdflatex ohyou
	pdflatex ohyou

tidy:
	rm -f ohyou.{aux,bbl,blg,log}

clean: tidy
	rm -f ohyou.pdf

.PHONY: all tidy clean
