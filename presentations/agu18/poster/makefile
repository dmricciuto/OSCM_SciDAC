SHELL=/bin/sh

main     = poster.tex
bibfiles = 
subfiles = 

default: doc
all: doc

.IGNORE:

clean:
	/bin/rm -rf $(main:.tex=.dvi) $(main:.tex=.log) $(main:.tex=.aux)  
	/bin/rm -rf $(main:.tex=.lof) $(main:.tex=.lot) $(main:.tex=.toc) 
	/bin/rm -rf $(main:.tex=.pdf) $(main:.tex=.ps)  $(main:.tex=.out)
	/bin/rm -rf $(main:.tex=.bbl) $(main:.tex=.blg) $(main:.tex=.brf)
	/bin/rm -rf core *~
	/bin/rm -rf *fls *fdb_latexmk *synctex.gz

doc: $(main) $(subfiles) $(bibfiles)
	pdflatex $(main)

# docB: $(main) $(subfiles) $(bibfiles)
# 	pdflatex $(main)
# 	bibtex   $(main:.tex=)
# 	pdflatex $(main)
# 	pdflatex $(main)

