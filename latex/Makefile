all: *.tex
#	pdflatex ps-ass1Rep.tex 
#	pdflatex ps-ass1Rep.tex 
#	pdflatex ps-ass1Rep.tex 
	latex ps-ass1Rep.tex
	bibtex ps-ass1Rep
	latex ps-ass1Rep.tex
	latex ps-ass1Rep.tex
	dvipdf ps-ass1Rep.dvi
	cp ps-ass1Rep.pdf ../
	xpdf ps-ass1Rep.pdf

clean:
	rm -rf *.toc *.aux *.log
	
