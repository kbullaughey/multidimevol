DOC = bullaughey-2012_evolution
SUP = $(DOC)-supporting_information

all: $(DOC).pdf $(SUP).pdf

$(DOC).pdf: $(DOC).tex $(DOC).bib mychicago.bst
	latex $(DOC)
	bibtex $(DOC)
	latex $(DOC)
	latex $(DOC)
	dvips -t letter -o $(DOC).ps $(DOC).dvi
	gs -sDEVICE=pdfwrite -sOutputFile=$(DOC).pdf -dNOPAUSE -dBATCH -sPAPERSIZE=letter -dAutoRotatePages=/None $(DOC).ps

$(SUP).pdf: $(SUP).tex
	pdflatex $(SUP)
	pdflatex $(SUP)

clean: 
	-rm $(DOC).blg $(DOC).bbl $(DOC).aux $(DOC).log $(DOC).pdf $(SUP).aux $(SUP).log $(DOC).ps $(DOC).dvi $(SUP).pdf
