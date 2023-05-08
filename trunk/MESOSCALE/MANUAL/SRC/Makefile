PROJ = user_manual

dvi: $(PROJ).tex *.tex
	latex -interaction nonstopmode $(PROJ).tex

all_format: dvi pdf ps

clean:
	-rm -f $(PROJ).log $(PROJ).aux $(PROJ).dvi $(PROJ).mtc* $(PROJ).bbl $(PROJ).blg

run: $(PROJ).dvi
	xdvi $(PROJ).dvi &

#pdf: $(PROJ).dvi
#	dvipdf $(PROJ).dvi

check: $(PROJ).tex *.tex
	make clean
	pdflatex $(PROJ).tex
	pdflatex $(PROJ).tex

pdf: $(PROJ).tex *.tex
	make clean
	pdflatex $(PROJ).tex
	bibtex $(PROJ)
	pdflatex $(PROJ).tex
	pdflatex $(PROJ).tex
	#xpdf $(PROJ).pdf

ps: $(PROJ).dvi
	dvips -t letter $(PROJ).dvi

#all: $(PROJ).tex
#	cp /u/forget/TXT/BIBLIO/newfred.bib newfred.bib_backup
#	pdflatex $(PROJ)
#	bibtex $(PROJ)
#	pdflatex $(PROJ)
#	pdflatex $(PROJ)
#	cp $(PROJ).pdf /u/aslmd/WWW/antichambre
#	#acroread $(PROJ).pdf
