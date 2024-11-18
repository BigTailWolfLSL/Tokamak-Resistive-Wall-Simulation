# General PDF compilation

all: MPhil_Poster_Template_A1.pdf MPhil_Poster_Template.pdf CDT_Poster_Template_A1.pdf CDT_Poster_Template.pdf

%.pdf : %.tex
	latexmk -pdf -pdflatex="pdflatex -interaction=nonstopmode" -use-make $<

clean:
	rm -f *.pdf *.aux *.bbl *.blg *.log *.nav *.out *.run.xml *.snm *.toc
