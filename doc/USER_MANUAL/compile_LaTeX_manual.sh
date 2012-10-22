#!/bin/sh

/bin/rm -rf *.dvi >  /dev/null
/bin/rm -rf *.log >  /dev/null
/bin/rm -rf *.out >  /dev/null
/bin/rm -rf *.aux >  /dev/null
/bin/rm -rf *.toc >  /dev/null
/bin/rm -rf *.blg >  /dev/null
/bin/rm -rf *.bbl >  /dev/null
/bin/rm -rf *.lof >  /dev/null
/bin/rm -rf *.lot >  /dev/null
/bin/rm -rf *.plt >  /dev/null
/bin/rm -rf *.fff >  /dev/null
/bin/rm -rf *.ttt >  /dev/null
/bin/rm -rf *.tit >  /dev/null
/bin/rm -rf *.spl >  /dev/null

pdflatex manual_SPECFEM3D_Cartesian
bibtex manual_SPECFEM3D_Cartesian
pdflatex manual_SPECFEM3D_Cartesian
pdflatex manual_SPECFEM3D_Cartesian
pdflatex manual_SPECFEM3D_Cartesian
pdflatex manual_SPECFEM3D_Cartesian

/bin/rm -rf *.dvi >  /dev/null
/bin/rm -rf *.log >  /dev/null
/bin/rm -rf *.out >  /dev/null
/bin/rm -rf *.aux >  /dev/null
/bin/rm -rf *.toc >  /dev/null
/bin/rm -rf *.blg >  /dev/null
/bin/rm -rf *.bbl >  /dev/null
/bin/rm -rf *.lof >  /dev/null
/bin/rm -rf *.lot >  /dev/null
/bin/rm -rf *.plt >  /dev/null
/bin/rm -rf *.fff >  /dev/null
/bin/rm -rf *.ttt >  /dev/null
/bin/rm -rf *.tit >  /dev/null
/bin/rm -rf *.spl >  /dev/null

