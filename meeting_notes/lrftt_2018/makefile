SHELL = /bin/csh
target = note
main = $(target).tex
##subfiles = *.tex bio/*.tex *.bib fig/*

all: pdf
default: pdf

.IGNORE:

dvi: $(target).dvi

ps: $(target).ps
pdf: $(target).pdf

clean:
	/bin/rm -f $(target).dvi $(target).log $(target).ps
	/bin/rm -f $(target).bbl $(target).aux $(target).blg
	/bin/rm -f $(target).lof $(target).fff $(target).toc
	/bin/rm -f $(target).pdf $(target).synctex.gz core

$(target).pdf: $(main) $(subfiles)
	pdflatex $(main)
#	bibtex   $(main:.tex=)
#	pdflatex $(main)
#	pdflatex $(main)
