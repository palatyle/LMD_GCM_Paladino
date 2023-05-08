#!/bin/bash
#  This script simly run pdflatex to make the documentation

# 0. cleanup from previous compilations
\rm -f main.pdf main.bbl main.aux main.blg main.idx main.log main.toc

# 1. make some files "on the fly" based on current deftank examples:
cd input
echo '\begin{verbatim}' > z2sig.tex
cat ../../deftank/z2sig.def >> z2sig.tex
echo '\end{verbatim}' >> z2sig.tex

echo '\begin{verbatim}' > run.tex
cat ../../deftank/run.def.64x48x49.MCD5 >> run.tex
echo '\end{verbatim}' >> run.tex

echo '\begin{verbatim}' > callphys.tex
cat ../../deftank/callphys.def.MCD5 >> callphys.tex
echo '\end{verbatim}' >> callphys.tex

echo '\begin{verbatim}' > run.def.1d.tex
cat ../../deftank/run.def.1d >> run.def.1d.tex
echo '\end{verbatim}' >> run.def.1d.tex

# 2. proceed with pdflatex
cd -
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex

# 3. copy output main.pdf to ../user_manual.pdf
\mv -f main.pdf ../user_manual.pdf
 
