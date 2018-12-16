
SHELL = /bin/csh

master = slides.tex 
target = $(master:.tex=)

files = *.tex

default: pdf
all: pdf

.IGNORE:

pdf: $(master:.tex=.pdf)

clean:
	/bin/rm -f $(target).aux
	/bin/rm -f $(target).log
	/bin/rm -f $(target).nav
	/bin/rm -f $(target).out
	/bin/rm -f $(target).pdf
	/bin/rm -f $(target).snm
	/bin/rm -f $(target).toc
	/bin/rm -f $(target).dvi
	/bin/rm -f core

$(master:.tex=.pdf): $(files)
	pdflatex $(master)
	pdflatex $(master)

