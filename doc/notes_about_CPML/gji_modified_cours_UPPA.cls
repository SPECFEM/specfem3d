% Documentclass gji.cls
% for Geophysical Journal International
% Based on the MN journal document style  
% 

% Version 2.0 for Latex 2e extensively rewritten by
% Brian Kennett, Pacific Region Office, GJI
 
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{gji_modified_cours_UPPA}[1995/11/30 v2.0 LaTeX2e document class]

\def\f@s@s{}% 
% 
 

% 
\typeout{Document Class `GJI' \f@s@s V2.0, released 1998 Nov} 
 
\def\f@ce{cm} 

%  Journals use two-sided printing. 
% 
\@twosidetrue 
\@twocolumntrue 
\@mparswitchtrue 
\def\ds@draft{\overfullrule 5pt} 
\def\ds@twocolumn{\@twocolumntrue} 
\def\ds@onecolumn{\@twocolumnfalse} 
 
% fasttrack option 
% 
\newif\if@fasttrack 
\def\ds@fasttrack{\@fasttracktrue} 
 
% landscape option 
% 
\newif\ifGJ@landscape 
\def\ds@landscape{\GJ@landscapetrue} 
 
% galley option 
% 
\newif\ifGJ@galley 
\def\ds@galley{\GJ@galleytrue} 
 
% referee option 
% 
\newif\ifGJ@referee 
\def\ds@referee{% 
 \GJ@refereetrue 
 \@twocolumnfalse 
} 
\newif\ifGJ@mreferee 
\def\ds@mreferee{% 
 \GJ@mrefereetrue 
 \@twocolumnfalse 
} 
\newif\ifGJ@camera 
\def\ds@camera{% 
 \GJ@cameratrue 
} 
\newif\ifGJ@extra 
\def\ds@extra{% 
 \GJ@extratrue 
} 
\newif\ifdoublesp@cing 
\def\ds@doublespacing{% 
 \doublesp@cingtrue 
}% 
\ExecuteOptions{a4paper,twoside,final}
\ProcessOptions
% 
%  **************************************** 
%  *               FONTS                  * 
%  **************************************** 
% 
\lineskip 1pt 
\normallineskip 1pt 
\ifGJ@referee 
  \renewcommand\baselinestretch{1.75}
\else
  \ifdoublesp@cing 
   \renewcommand\baselinestretch{2} 
 \else 
   \renewcommand\baselinestretch{1} 
\fi\fi 
\ifGJ@mreferee 
  \renewcommand\baselinestretch{1.75}
\fi

\renewcommand\normalsize{\@setfontsize\normalsize\@ixpt{12}%
   \abovedisplayskip 4pt plus 1pt minus 1pt 
   \belowdisplayskip \abovedisplayskip 
   \abovedisplayshortskip 4pt plus 1pt 
   \belowdisplayshortskip \abovedisplayshortskip 
   \let\@listi\@listI} 
\newcommand\small{\@setfontsize\small\@viiipt{10}% 
 \abovedisplayskip 4pt plus 1pt minus 1pt 
 \belowdisplayskip \abovedisplayskip 
 \abovedisplayshortskip 4pt plus 1pt 
 \belowdisplayshortskip \abovedisplayshortskip 
 \def\@listi{\leftmargin\leftmargini 
  \topsep 2pt plus 1pt minus 1pt 
  \parsep \z@ 
  \itemsep 2pt}} 
\newcommand\footnotesize{\@setfontsize\footnotesize\@viiipt{10}% 
 \abovedisplayskip 4pt plus 1pt minus 1pt 
 \belowdisplayskip \abovedisplayskip 
 \abovedisplayshortskip 4pt plus 1pt 
 \belowdisplayshortskip \abovedisplayshortskip 
 \renewcommand\@listi{\leftmargin\leftmargini 
  \topsep 2pt plus 1pt minus 1pt 
  \parsep \z@ 
  \itemsep 2pt}} 
\newcommand\scriptsize{\@setfontsize\scriptsize\@vipt{8}} 
\newcommand\tiny{\@setfontsize\tiny\@vpt{6}} 

\newcommand\large{\@setfontsize\large\@xpt{12}} 
\newcommand\Large{\@setfontsize\Large\@xiipt{14}} 
\newcommand\LARGE{\@setfontsize\LARGE\@xivpt{17}} 
\newcommand\huge{\@setfontsize\huge\@xviipt{20}} 
\newcommand\Huge{\@setfontsize\Huge\@xxpt{25}} 

\if@twocolumn 
   \renewcommand\normalsize{\@setfontsize\normalsize\@ixpt{11}%
   \abovedisplayskip 6pt plus 2pt minus 2pt 
   \belowdisplayskip \abovedisplayskip 
   \abovedisplayshortskip 6pt plus 2pt 
   \belowdisplayshortskip \abovedisplayshortskip 
   \let\@listi\@listI}
  \renewcommand\large{\@setfontsize\large\@xpt{11}}
\fi 
 
\ifGJ@referee 
  \renewcommand\normalsize{\@setfontsize\normalsize\@xiipt{14}%
   \abovedisplayskip 8pt plus 1pt minus 1pt 
   \belowdisplayskip \abovedisplayskip 
   \abovedisplayshortskip 4pt plus 1pt 
   \belowdisplayshortskip \abovedisplayshortskip 
   \let\@listi\@listI}
 \renewcommand\small{\@setfontsize\small\@xipt{11}% 
  \abovedisplayskip 5pt plus 1pt minus 1pt 
  \belowdisplayskip \abovedisplayskip 
  \abovedisplayshortskip 4pt plus 1pt 
  \belowdisplayshortskip \abovedisplayshortskip 
 \def\@listi{\leftmargin\leftmargini 
  \topsep 2pt plus 1pt minus 1pt 
  \parsep \z@ 
  \itemsep 2pt}}  
  \renewcommand\large{\@setfontsize\large\@xiipt{13}} 
  \renewcommand\Large{\@setfontsize\Large\@xiipt{15}}
  \renewcommand\LARGE{\@setfontsize\LARGE\@xivpt{18}}
\fi
 
\ifGJ@mreferee 
  \renewcommand\normalsize{\@setfontsize\normalsize\@xipt{11}%
   \abovedisplayskip 8pt plus 1pt minus 1pt 
   \belowdisplayskip \abovedisplayskip 
   \abovedisplayshortskip 4pt plus 1pt 
   \belowdisplayshortskip \abovedisplayshortskip 
   \let\@listi\@listI} 
 \renewcommand\small{\@setfontsize\small\@xpt{10}% 
  \abovedisplayskip 5pt plus 1pt minus 1pt 
  \belowdisplayskip \abovedisplayskip 
  \abovedisplayshortskip 4pt plus 1pt 
  \belowdisplayshortskip \abovedisplayshortskip 
 \def\@listi{\leftmargin\leftmargini 
  \topsep 2pt plus 1pt minus 1pt 
  \parsep \z@ 
  \itemsep 2pt}} 
  \renewcommand\large{\@setfontsize\large\@xiipt{13}} 
  \renewcommand\Large{\@setfontsize\Large\@xiipt{15}}
  \renewcommand\LARGE{\@setfontsize\LARGE\@xivpt{18}}
\fi 

\ifGJ@camera
   \renewcommand\normalsize{\@setfontsize\normalsize\@xpt{11}% 
   \abovedisplayskip 6pt plus 2pt minus 2pt 
   \belowdisplayskip \abovedisplayskip 
   \abovedisplayshortskip 6pt plus 2pt 
   \belowdisplayshortskip \abovedisplayshortskip 
   \let\@listi\@listI} 
   \renewcommand\large{\@setfontsize\large\@xipt{12}}
\fi

\normalsize 
 
\newcommand{\mitbf}[1]{  
  \hbox{\mathversion{bold}$#1$}}  
  \newcommand{\rmn}[1] {{\mathrm #1}}  
  \newcommand{\itl}[1] {{\mathit #1}}  
  \newcommand{\bld}[1] {{\mathbf #1}}  
 
\DeclareOldFontCommand{\rm}{\normalfont\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\normalfont\sffamily}{\mathsf}
\DeclareOldFontCommand{\tt}{\normalfont\ttfamily}{\mathtt}
\DeclareOldFontCommand{\bf}{\normalfont\bfseries}{\mathbf}
\DeclareOldFontCommand{\it}{\normalfont\itshape}{\mathit}
\DeclareOldFontCommand{\sl}{\normalfont\slshape}{\@nomath\sl}
\DeclareOldFontCommand{\sc}{\normalfont\scshape}{\@nomath\sc}
\DeclareRobustCommand\cal{\@fontswitch{\relax}{\mathcal}}
\DeclareRobustCommand\mit{\@fontswitch{\relax}{\mathchoice}}
%\DeclareRobustCommand*\cal{\@fontswitch{\relax}{\mathcal}}
%\DeclareRobustCommand*\mit{\@fontswitch{\relax}{\mathchoice}}
%

 
% Astronomy and Astrophysics symbol macros 
 
\newcommand\la{\mathrel{\mathchoice {\vcenter{\offinterlineskip\halign{\hfil 
$\reset@font\displaystyle##$\hfil\cr<\cr\sim\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\textstyle##$\hfil\cr 
<\cr\sim\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\scriptstyle##$\hfil\cr 
<\cr\sim\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\scriptscriptstyle##$\hfil\cr 
<\cr\sim\cr}}}}} 
 
\newcommand\ga{\mathrel{\mathchoice {\vcenter{\offinterlineskip\halign{\hfil 
$\reset@font\displaystyle##$\hfil\cr>\cr\sim\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\textstyle##$\hfil\cr 
>\cr\sim\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\scriptstyle##$\hfil\cr 
>\cr\sim\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\scriptscriptstyle##$\hfil\cr 
>\cr\sim\cr}}}}} 
 
\newcommand\getsto{\mathrel{\mathchoice {\vcenter{\offinterlineskip 
\halign{\hfil 
$\reset@font\displaystyle##$\hfil\cr\gets\cr\to\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\textstyle##$\hfil\cr\gets 
\cr\to\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\scriptstyle##$\hfil\cr\gets 
\cr\to\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\scriptscriptstyle##$\hfil\cr 
\gets\cr\to\cr}}}}} 
 
\newcommand\cor{\mathrel{\mathchoice {\hbox{$\widehat=$}}{\hbox{$\widehat=$}} 
{\hbox{$\reset@font\scriptstyle\hat=$}} 
{\hbox{$\reset@font\scriptscriptstyle\hat=$}}}} 
 
\newcommand\lid{\mathrel{\mathchoice {\vcenter{\offinterlineskip\halign{\hfil 
$\reset@font\displaystyle##$\hfil\cr<\cr\noalign{\vskip1.2pt}=\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\textstyle##$\hfil\cr<\cr 
\noalign{\vskip1.2pt}=\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\scriptstyle##$\hfil\cr<\cr 
\noalign{\vskip1pt}=\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\scriptscriptstyle##$\hfil\cr 
<\cr 
\noalign{\vskip0.9pt}=\cr}}}}} 
 
\newcommand\gid{\mathrel{\mathchoice {\vcenter{\offinterlineskip\halign{\hfil 
$\reset@font\displaystyle##$\hfil\cr>\cr\noalign{\vskip1.2pt}=\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\textstyle##$\hfil\cr>\cr 
\noalign{\vskip1.2pt}=\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\scriptstyle##$\hfil\cr>\cr 
\noalign{\vskip1pt}=\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\scriptscriptstyle##$\hfil\cr 
>\cr 
\noalign{\vskip0.9pt}=\cr}}}}} 
 
\newcommand\sol{\mathrel{\mathchoice {\vcenter{\offinterlineskip\halign{\hfil 
$\reset@font\displaystyle##$\hfil\cr\sim\cr<\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\textstyle##$\hfil\cr\sim\cr 
<\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\scriptstyle##$\hfil\cr\sim\cr 
<\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\scriptscriptstyle##$\hfil\cr 
\sim\cr<\cr}}}}} 
 
\newcommand\sog{\mathrel{\mathchoice {\vcenter{\offinterlineskip\halign{\hfil 
$\reset@font\displaystyle##$\hfil\cr\sim\cr>\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\textstyle##$\hfil\cr\sim\cr 
>\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\scriptstyle##$\hfil\cr 
\sim\cr>\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\scriptscriptstyle##$\hfil\cr 
\sim\cr>\cr}}}}} 
 
\newcommand\lse{\mathrel{\mathchoice {\vcenter{\offinterlineskip\halign{\hfil 
$\reset@font\displaystyle##$\hfil\cr<\cr\simeq\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\textstyle##$\hfil\cr 
<\cr\simeq\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\scriptstyle##$\hfil\cr 
<\cr\simeq\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\scriptscriptstyle##$\hfil\cr 
<\cr\simeq\cr}}}}} 
 
\newcommand\gse{\mathrel{\mathchoice {\vcenter{\offinterlineskip\halign{\hfil 
$\reset@font\displaystyle##$\hfil\cr>\cr\simeq\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\textstyle##$\hfil\cr 
>\cr\simeq\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\scriptstyle##$\hfil\cr 
>\cr\simeq\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\scriptscriptstyle##$\hfil\cr 
>\cr\simeq\cr}}}}} 
 
\newcommand\grole{\mathrel{\mathchoice {\vcenter{\offinterlineskip\halign{\hfil 
$\reset@font\displaystyle##$\hfil\cr>\cr\noalign{\vskip-1.5pt}<\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\textstyle##$\hfil\cr 
>\cr\noalign{\vskip-1.5pt}<\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\scriptstyle##$\hfil\cr 
>\cr\noalign{\vskip-1pt}<\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\scriptscriptstyle##$\hfil\cr 
>\cr\noalign{\vskip-0.5pt}<\cr}}}}} 
 
\newcommand\leogr{\mathrel{\mathchoice {\vcenter{\offinterlineskip\halign{\hfil 
$\reset@font\displaystyle##$\hfil\cr<\cr\noalign{\vskip-1.5pt}>\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\textstyle##$\hfil\cr 
<\cr\noalign{\vskip-1.5pt}>\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\scriptstyle##$\hfil\cr 
<\cr\noalign{\vskip-1pt}>\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\scriptscriptstyle##$\hfil\cr 
<\cr\noalign{\vskip-0.5pt}>\cr}}}}} 
 
\newcommand\loa{\mathrel{\mathchoice {\vcenter{\offinterlineskip\halign{\hfil 
$\reset@font\displaystyle##$\hfil\cr<\cr\approx\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\textstyle##$\hfil\cr 
<\cr\approx\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\scriptstyle##$\hfil\cr 
<\cr\approx\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\scriptscriptstyle##$\hfil\cr 
<\cr\approx\cr}}}}} 
 
\newcommand\goa{\mathrel{\mathchoice {\vcenter{\offinterlineskip\halign{\hfil 
$\reset@font\displaystyle##$\hfil\cr>\cr\approx\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\textstyle##$\hfil\cr 
>\cr\approx\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\scriptstyle##$\hfil\cr 
>\cr\approx\cr}}} 
{\vcenter{\offinterlineskip\halign{\hfil$\reset@font\scriptscriptstyle##$\hfil\cr 
>\cr\approx\cr}}}}} 
 
\newcommand\sun{\hbox{$\odot$}} 
\newcommand\earth{\hbox{$\oplus$}}
\newcommand\degr{\hbox{$^\circ$}} 
\newcommand\diameter{{\ifmmode\mathchoice 
{\ooalign{\hfil\hbox{$\reset@font\displaystyle/$}\hfil\crcr 
{\hbox{$\reset@font\displaystyle\mathchar"20D$}}}} 
{\ooalign{\hfil\hbox{$\reset@font\textstyle/$}\hfil\crcr 
{\hbox{$\reset@font\textstyle\mathchar"20D$}}}} 
{\ooalign{\hfil\hbox{$\reset@font\scriptstyle/$}\hfil\crcr 
{\hbox{$\reset@font\scriptstyle\mathchar"20D$}}}} 
{\ooalign{\hfil\hbox{$\reset@font\scriptscriptstyle/$}\hfil\crcr 
{\hbox{$\reset@font\scriptscriptstyle\mathchar"20D$}}}} 
\else{\ooalign{\hfil/\hfil\crcr\mathhexbox20D}}% 
\fi}} 
 
\newcommand\sq{\ifmmode\squareforqed\else{\unskip\nobreak\hfil 
\penalty50\hskip1em\null\nobreak\hfil\squareforqed 
\parfillskip=0pt\finalhyphendemerits=0\endgraf}\fi} 
\newcommand\squareforqed{\hbox{\rlap{$\sqcap$}$\sqcup$}} 
 
\newcommand{\romn}[1] {{\mathrm #1}} 
 
\newcommand\fd{\hbox{$.\!\!^{\reset@font\romn d}$}} 
\newcommand\fh{\hbox{$.\!\!^{\reset@font\romn h}$}} 
\newcommand\fm{\hbox{$.\!\!^{\reset@font\romn m}$}} 
\newcommand\fs{\hbox{$.\!\!^{\reset@font\romn s}$}} 
\newcommand\fdg{\hbox{$.\!\!^\circ$}} 
\newcommand\farcm{\hbox{$.\mkern-4mu^\prime$}} 
\newcommand\farcs{\hbox{$.\!\!^{\prime\prime}$}} 
\newcommand\fp{\hbox{$.\!\!^{\reset@font\reset@font\scriptscriptstyle\romn p}$}} 
\newcommand\arcmin{\hbox{$^\prime$}} 
\newcommand\arcsec{\hbox{$^{\prime\prime}$}} 

\def\micron{\hbox{$\mu$m}}

\newcommand\case[2]{\hbox{$\frac{#1}{#2}$}}
\newcommand\slantfrac[2]{\hbox{$\,^#1\!/_#2$}}
\newcommand\onehalf{\slantfrac{1}{2}}
\newcommand\onethird{\slantfrac{1}{3}}
\newcommand\twothirds{\slantfrac{2}{3}}
\newcommand\onequarter{\slantfrac{1}{4}}
\newcommand\threequarters{\slantfrac{3}{4}}
\newcommand\fhf{\case{1}{2}} 
\newcommand\fqr{\case{1}{4}}
% 
%\chemical from Kopka & Daly
\newlength\fntxvi 
\newlength\fntxvii
\newcommand{\chemical}[1]
{{\fontencoding{OMS}\fontfamily{cmsy}\selectfont
  \fntxvi\the\fontdimen16\font
  \fntxvii\the\fontdimen17\font
  \fontdimen16\font=2pt\fontdimen17\font=2pt
  $\mathrm{#1}$
  \fontencoding{OMS}\fontfamily{cmsy}\selectfont
  \fontdimen16\font=\fntxvi \fontdimen17\font=\fntxvii}}
 
%  **************************************** 
%  *            PAGE LAYOUT               * 
%  **************************************** 
% 
% SIDE MARGINS: 
\if@twocolumn 
  \setlength\oddsidemargin{-2.5pc}
  \setlength\evensidemargin{-2.5pc}
\else 
 \ifGJ@referee 
  \setlength\oddsidemargin{\z@}
  \setlength\evensidemargin{\z@} 
 \else 
  \setlength\oddsidemargin{-2pc}
  \setlength\evensidemargin{-2pc} 
 \fi 
\fi 
\setlength\marginparwidth{4pc}
\setlength\marginparsep{0.5pc}
 
% VERTICAL SPACING: 
%\if@twocolumn 
  \setlength\topmargin{-2pc} 
%\else 
  \setlength\topmargin{\z@}  
%\fi 
\setlength{\headheight}{14pt} 
\setlength{\headsep}   {15pt} 
\setlength{\topskip}   {9pt}
%\setlength{\footheight}{11pt} 
\setlength{\footskip}  {26pt} 
% 
% DIMENSION OF TEXT:  
% 
\newdimen\GJ@measure 
\GJ@measure 42pc 
% 
\setlength\textheight{682pt}
\setlength\textwidth{\GJ@measure} 
% 
\ifGJ@landscape 
 \setlength\textwidth{\textheight} 
 \setlength\textheight{\GJ@measure}
\fi 
\ifGJ@referee 
 \setlength\textwidth{39pc} 
\fi 
\setlength\columnsep{2pc} 
\setlength\columnseprule{0pt}

\ifGJ@mreferee 
\setlength {\textwidth}     {156.0mm}
\setlength {\textheight}    {230.0mm}
\setlength {\evensidemargin}{0mm}
\setlength {\oddsidemargin} {0mm}
\setlength {\topmargin} {0mm}
\setlength {\baselineskip} {17pt}
\setlength {\parindent}     {5mm}
\fi
 
% FOOTNOTES: 
\setlength\footnotesep{6pt} 
\setlength\skip{\footins 19.5pt plus 12pt minus 1pt} 
 
% FLOATS: 
\setlength\floatsep{ 12pt plus 5.5pt minus 1pt} 
\setlength\textfloatsep{ 12pt plus 7.5pt}  
\setlength\intextsep{ 12pt plus 5.5pt minus 2pt} 
%%\setlength\@maxsep{ 15pt }
\setlength\dblfloatsep{ 12pt plus 5.5pt minus 2pt }
\setlength\dbltextfloatsep{ 15pt plus 4.5pt minus 3pt} 
%%\setlength\@dblmaxsep{ 15pt }
\setlength\@fptop{ 0pt plus 0fil} 
\setlength\@fpsep{ 12pt plus 0fil }
\setlength\@fpbot{ 0pt plus 3fil }
\setlength\@dblfptop{ 0pt plus 0fil} 
\setlength\@dblfpsep{ 12pt plus 0fil} 
\setlength\@dblfpbot{ 0pt plus 3fil} 
\setlength\marginparpush{ 5pt }
 
%  **************************************** 
%  *           PARAGRAPHING               * 
%  **************************************** 
% 
\setlength\parskip{ \z@ plus .1pt }
\setlength\parindent{ 18pt }
\setlength\partopsep{ 0pt plus 1pt }
\@lowpenalty   51 
\@medpenalty  151 
\@highpenalty 301 
\@beginparpenalty -\@lowpenalty 
\@endparpenalty -\@lowpenalty 
\@itempenalty -\@lowpenalty 
\clubpenalty=0 
\widowpenalty=10000 
 
%  **************************************** 
%  *       PARTS & SECTIONS NUMBERS       * 
%  **************************************** 
% 
\newcounter {part} 
\newcounter {section} 
\newcounter {subsection}[section] 
\newcounter {subsubsection}[subsection] 
\newcounter {paragraph}[subsubsection] 
\newcounter {subparagraph}[paragraph] 
\renewcommand\thepart {\arabic{part}} 
\renewcommand\thesection {\arabic{section}} 
\renewcommand\thesubsection    {\thesection.\arabic{subsection}} 
\renewcommand\thesubsubsection {\thesubsection.\arabic{subsubsection}} 
\renewcommand\theparagraph     {\thesubsubsection.\arabic{paragraph}} 
\renewcommand\thesubparagraph  {\theparagraph.\arabic{subparagraph}} 

% Sometimes people want to number equations by section, rather than
% sequentially through the whole paper.

\newcommand\eqsecnum{
    \@newctr{equation}[section]
   \def\theequation{\normalsize\mbox{\arabic{section}.\arabic{equation}}}}   
%    \def\theequation{\mbox{\arabic{section}.\arabic{equation}}}}
 
\newcommand\eqsubsecnum{
    \@newctr{equation}[subsection]
    \def\theequation{\normalsize%
    \mbox{\arabic{section}.\arabic{subsection}.\arabic{equation}}}}
%    \def\theequation
%   {\mbox{\arabic{section}.\arabic{subsection}.\arabic{equation}}}}
 

%  **************************************** 
%  *                PARTS                 * 
%  **************************************** 
% 
\def\part{\par \addvspace{4ex}\@afterindentfalse 
 \secdef\@part\@spart} 
\def\@part[#1]#2{\ifnum \c@secnumdepth >\m@ne 
  \refstepcounter{part} 
  \addcontentsline{toc}{part}{Part \thepart: #1} 
 \else \addcontentsline{toc}{part}{#1} 
 \fi 
 {\parindent 0pt \raggedright 
  \ifnum \c@secnumdepth >\m@ne 
   \reset@font\large\rm PART 
   \ifcase\thepart \or ONE \or TWO \or THREE \or FOUR \or FIVE 
    \or SIX \or SEVEN \or EIGHT \or NINE \or TEN \else \fi 
   \par \nobreak 
  \fi 
  \reset@font\LARGE \rm #2 \markboth{}{}\par } 
 \nobreak \vskip 3ex \@afterheading} 
\def\@spart#1{{\parindent 0pt \raggedright 
  \reset@font\LARGE \rm #1\par} 
 \nobreak \vskip 3ex \@afterheading} 
 
%  **************************************** 
%  *               SECTIONS               * 
%  **************************************** 
% 
% \@startsection modified to leave less space before \subsection 
%    if it immediately follows \section 
 
% \@startsection {NAME}{LEVEL}{INDENT}{BEFORESKIP}{AFTERSKIP}{STYLE}  
 
\def\@startsection#1#2#3#4#5#6{\if@noskipsec \leavevmode \fi 
   \par \@tempskipa #4\relax 
   \@afterindenttrue 
   \ifdim \@tempskipa <\z@ \@tempskipa -\@tempskipa \@afterindentfalse\fi 
   \if@nobreak 
     \everypar{} \ifnum#2=2 \vskip 0pt plus1pt\fi % was 6pt 
   \else 
     \addpenalty{\@secpenalty}\addvspace{\@tempskipa} 
   \fi \@ifstar 
     {\@ssect{#3}{#4}{#5}{#6}}{\@dblarg{\@sect{#1}{#2}{#3}{#4}{#5}{#6}}}} 
 
%\def\GJ@hangraggedright{\rightskip\@flushglue \let\\=\newline} 
%\def\GJ@raggedright{\rightskip\@flushglue \let\\=\@centercr 
% \parindent=\GJ@indent} 
 
\def\GJ@hangraggedright{\rightskip=\@flushglue \let\\=\@centercr \parindent=0pt} 
 
\def\section{\@startsection{section}{1}{\z@} 
 {-24pt plus -12pt minus -1pt}{6pt} 
 {\GJ@hangraggedright\reset@font\normalsize\bf}} 
\def\subsection{\@startsection{subsection}{2}{\z@} 
 {-18pt plus -9pt minus -1pt}{6pt} 
 {\GJ@hangraggedright\reset@font\normalsize\bf}} 
\def\subsubsection{\@startsection{subsubsection}{3}{\z@} 
 {-18pt plus -9pt minus -1pt}{6pt} 
 {\GJ@hangraggedright\reset@font\normalsize\it}} 
\def\paragraph{\@startsection{paragraph}{4}{\z@} 
 {12pt plus 2.25pt minus 1pt}{-0.5em}{\reset@font\normalsize\bf}} 
\def\subparagraph{\@startsection{subparagraph}{5}{\parindent} 
 {12pt plus 2.25pt minus 1pt}{-0.5em}{\reset@font\normalsize\it}} 
\setcounter{secnumdepth}{4} 
% 
\newif\if@levelone 
% 
% \@sect{NAME}{LEVEL}{INDENT}{BEFORESKIP}{AFTERSKIP}{STYLE}[ARG1]{ARG2}  
\def\@sect#1#2#3#4#5#6[#7]#8{% 
  \ifnum #2>\c@secnumdepth 
  \def\@svsec{}% 
 \else 
  \refstepcounter{#1} 
   \ifnum #2=\@ne \global\@levelonetrue % if level=1 ie. section 
   \ifGJ@appendix \edef\@svsec{}% 
     \else \edef\@svsec{\csname the#1\endcsname\hskip 1em}% 
   \fi 
  \else \edef\@svsec{\csname the#1\endcsname\hskip 1em}% 
  \fi 
 \fi 
 \@tempskipa #5\relax 
 \ifdim \@tempskipa>\z@ 
  \begingroup #6\relax 
   \ifnum #2=\@ne 
    \ifGJ@appendix 
     \@hangfrom{\hskip #3\relax\@svsec}{\interlinepenalty \@M  
      APPENDIX \csname the#1\endcsname:\hskip 0.5em\uppercase{#8}\par}% 
    \else 
     \@hangfrom{\hskip #3\relax\@svsec}{\interlinepenalty \@M  
      \uppercase{#8}\par}% 
    \fi 
   \else 
    \@hangfrom{\hskip #3\relax\@svsec}{\interlinepenalty \@M #8\par}% 
   \fi 
  \endgroup 
   \csname #1mark\endcsname{#7}% 
   \addcontentsline{toc}{#1}{\ifnum #2>\c@secnumdepth \else %ajw 
   \protect\numberline{\csname the#1\endcsname}\fi #7}%ajw 
 \else 
  \def\@svsechd{#6\hskip #3\@svsec \ifnum #2=\@ne\uppercase{#8}\else #8\fi 
  \csname #1mark\endcsname{#7} 
  \addcontentsline{toc}{#1}{\ifnum #2>\c@secnumdepth \else %ajw 
  \protect\numberline{\csname the#1\endcsname}\fi#7}%ajw 
  }\fi 
 \@xsect{#5}} 
 
% from latex.tex 
%   \@ssect{INDENT}{BEFORESKIP}{AFTERSKIP}{STYLE}{ARG}  
\def\@ssect#1#2#3#4#5{% 
  \@tempskipa #3\relax 
    \ifdim \@tempskipa>\z@ 
      \begingroup #4\@hangfrom{\hskip #1}{% 
      \interlinepenalty \@M  
      \if@levelone 
        \uppercase{#5}\global\@levelonefalse 
      \else{#5} 
      \fi\par}% 
     \endgroup 
    \else  
      \def\@svsechd{#4\hskip #1\relax 
      \if@levelone 
        \uppercase{#5}\global\@levelonefalse 
      \else{#5} 
      \fi} 
    \fi 
  \@xsect{#3}% 
} 
 
\newif\ifGJ@appendix 
\def\appendix{\par 
 \GJ@appendixtrue 
 \setcounter{section}{0} 
% \def\thesection{A\arabic{section}} 
 \renewcommand\thesection{\Alph{section}} 
 \renewcommand\thesubsection{\thesection\arabic{subsection}} 
 \setcounter{equation}{0} 
% \def\theequation{A\arabic{equation}} 
 \@newctr{equation}[section]
\renewcommand\theequation{\hbox{\normalsize\Alph{section}.\arabic{equation}}}
%\renewcommand\theequation{\hbox{\Alph{section}.\arabic{equation}}}
\setcounter{figure}{0} 
 \renewcommand\thefigure{A\@arabic\c@figure} 
 \setcounter{table}{0} 
 \renewcommand\thetable{A\@arabic\c@table} 
} 

\newcommand\resetfigno{
 \setcounter{figure}{0} 
 \renewcommand\thefigure{\@arabic\c@figure} 
}
%  **************************************** 
%  *               SKIPS                  * 
%  **************************************** 
% 
\newskip\@indentskip 
\newskip\smallindent 
\newskip\@footindent 
\newskip\@leftskip 
\@indentskip=18pt 
\smallindent=1.5em 
\@footindent=\smallindent 
\@leftskip=\z@ 
 
%  **************************************** 
%  *                LISTS                 * 
%  **************************************** 
% 
\setlength\leftmargini{   \z@ }
\setlength\leftmarginii{  1.25pc }
\setlength\leftmarginiii{ 1.25pc }
\setlength\leftmarginiv{  1.25pc } 
\setlength\leftmarginv{  1pc }
\setlength\leftmarginvi{  1pc }
\setlength\leftmargin{\leftmargini} 
\setlength\labelsep{ 0.5em} 
\labelwidth\leftmargini\advance\labelwidth-\labelsep 
% 
% label macros for Range-Left and Range-Right labels 
\def\makeRLlabel#1{\rlap{#1}\hss} 
\def\makeRRlabel#1{\hss\llap{#1}} 
\def\makenewlabel#1{#1} 
% 
\def\@listI{\leftmargin\leftmargini 
 \parsep \z@ 
 \topsep 6pt plus 1pt minus 1pt 
 \itemsep \z@ plus .1pt 
} 
\let\@listi\@listI 
\@listi 
\def\@listii{\leftmargin\leftmarginii 
 \labelwidth\leftmarginii\advance\labelwidth-\labelsep 
 \topsep 6pt plus 1pt minus 1pt 
 \parsep \z@ 
 \itemsep \z@ plus .1pt 
} 
\def\@listiii{\leftmargin\leftmarginiii 
 \labelwidth\leftmarginiii\advance\labelwidth-\labelsep 
 \topsep 6pt plus 1pt minus 1pt 
 \parsep \z@ 
 \partopsep \z@ 
 \itemsep \topsep 
} 
\def\@listiv{\leftmargin\leftmarginiv 
 \labelwidth\leftmarginiv\advance\labelwidth-\labelsep 
} 
\def\@listv{\leftmargin\leftmarginv 
 \labelwidth\leftmarginv\advance\labelwidth-\labelsep 
} 
\def\@listvi{\leftmargin\leftmarginvi 
 \labelwidth\leftmarginvi\advance\labelwidth-\labelsep 
} 
% 
% ITEMIZE 
\def\itemize{% 
  \ifnum \@itemdepth >3 \@toodeep 
  \else \advance\@itemdepth \@ne 
   \edef\@itemitem{labelitem\romannumeral\the\@itemdepth}% 
   \list{\csname\@itemitem\endcsname}% 
  {\let\makelabel\makeRRlabel}% 
   \leftmargin 3.5em\labelsep 0.5em\itemindent 2.0em\labelwidth 0.5em% 
  \fi} 
\def\labelitemi{$\bullet$} 
\def\labelitemii{--} 
\def\labelitemiii{$\cdot$} 
\def\labelitemiv{*} 
% 
% ENUMERATE 
%  -- enumerate with optional argument to set left margin... 
\def\enumerate{\ifnum \@enumdepth >3 \@toodeep \else 
  \advance\@enumdepth \@ne 
  \edef\@enumctr{enum\romannumeral\the\@enumdepth}% 
 \fi 
 \@ifnextchar [{\@enumeratetwo}{\@enumerateone}% 
} 
\def\@enumeratetwo[#1]{% 
 \list{\csname label\@enumctr\endcsname}% 
  {\settowidth\labelwidth{[#1]} 
   \leftmargin\z@  
   \usecounter{\@enumctr} 
   \let\makelabel\makenewlabel} 
} 
\def\@enumerateone{% 
 \list{\csname label\@enumctr\endcsname}% 
  {\usecounter{\@enumctr}\itemindent=1em 
   \let\makelabel\makenewlabel}} 
% 
\renewcommand\theenumi{{\rm (\roman{enumi})}}
\renewcommand\theenumii{{\rm (\alph{enumii})}}
\renewcommand\theenumiii{{\rm \arabic{enumiii}.}}
\renewcommand\theenumiv{{\rm \Alph{enumiv}.}}

\newcommand\labelenumi{\theenumi}
\newcommand\labelenumii{\theenumii}
\newcommand\labelenumiii{\theenumiii}
\newcommand\labelenumiv{\theenumiv}

\renewcommand\p@enumii{\theenumi}
\renewcommand\p@enumiii{\theenumi(\theenumii)}
\renewcommand\p@enumiv{\p@enumiii\theenumiii}
% 
% DESCRIPTION 
%\def\descriptionlabel[1]{\hspace\labelsep\bfseries #1} 
%\def\description{\list{}{\labelwidth\z@ \itemindent-\leftmargin 
%  \leftmargin 1em 
%  \itemindent-1em 
%}} 
\newcommand\descriptionlabel[1]{\bfseries #1} 
\def\description{\list{}{\labelwidth\z@ 
  \leftmargin \z@ 
  \itemindent 1em 
}} 
\let\enddescription\endlist 
 
%  **************************************** 
%  *            DISPLAYED TEXT            * 
%  **************************************** 
% 
% VERSE 
\def\verse{\let\\=\@centercr 
 \list{}{\itemsep\z@ 
  \itemindent -\@indentskip 
  \listparindent \itemindent 
  \rightmargin\leftmargin 
  \advance\leftmargin \@indentskip}\item[]} 
\let\endverse\endlist 
% 
% QUOTATION 
\def\quotation{\list{}{\listparindent \smallindent 
%  \itemindent\listparindent 
  \leftmargin\z@\rightmargin\leftmargin 
  \parsep 0pt plus 1pt}\item[]\reset@font\small} 
\let\endquotation=\endlist 
% 
% QUOTE 
\def\quote{\list{}{\leftmargin\z@\rightmargin\leftmargin}\item[]\reset@font\small} 
\let\endquote=\endlist 
 
%  **************************************** 
%  *         OTHER ENVIRONMENTS           * 
%  **************************************** 
% 
% THEOREM 
\def\@begintheorem#1#2{\reset@font\rm \trivlist \item[\hskip \labelsep{\reset@font\bf #1\ #2.}]} 
\def\@opargbegintheorem#1#2#3{\reset@font\rm \trivlist 
  \item[\hskip \labelsep{\reset@font\bf #1\ #2.\ (#3)}]} 
\def\@endtheorem{\endtrivlist} 
% 
% PROOF 
\def\proof{\reset@font\rm \trivlist \item[\hskip \labelsep{\reset@font\it Proof.}]} 
\def\endproof{\hspace*{1em}{\begin{picture}(6.5,6.5)% 
  \put(0,0){\framebox(6.5,6.5){}}\end{picture}}\endtrivlist} 
\@namedef{proof*}{\reset@font\rm \trivlist \item[\hskip \labelsep{\reset@font\it Proof.}]} 
\@namedef{endproof*}{\endtrivlist} 
\def\proofbox{\begin{picture}(6.5,6.5)% 
  \put(0,0){\framebox(6.5,6.5){}}\end{picture}} 
% 
% TITLEPAGE 
\def\titlepage{\@restonecolfalse\if@twocolumn\@restonecoltrue\onecolumn 
  \else \newpage \fi \thispagestyle{empty}\c@page\z@} 
\def\endtitlepage{\if@restonecol\twocolumn \else \newpage \fi} 
 
% ARRAY AND TABULAR 
% 
\setlength\arraycolsep{   5pt } 
\setlength\tabcolsep{     6pt }
\setlength\arrayrulewidth{ .5pt }
\setlength\doublerulesep{ 1.5pt }
 
\def\tabular{\def\@halignto{} 
 \def\hline{\noalign{\ifnum0=`}\fi 
  \vskip 3pt 
  \hrule \@height \arrayrulewidth 
%   {\hbox to 0.667\textwidth{\hrulefill}} 
  \vskip 3pt 
  \futurelet \@tempa\@xhline} 
 \def\fullhline{\noalign{\ifnum0=`}\fi 
  \vskip 3pt 
  \hrule \@height \arrayrulewidth 
%   {\hbox to \textwidth{\hrulefill}} 
  \vskip 3pt 
  \futurelet \@tempa\@xhline} 
 \def\@xhline{\ifx\@tempa\hline 
   \vskip -6pt 
   \vskip \doublerulesep 
  \fi 
  \ifnum0=`{\fi}} 
  \def\@arrayrule{\@addtopreamble{\hskip -.5\arrayrulewidth 
%                                  \vrule \@width \arrayrulewidth 
                                  \hskip .5\arrayrulewidth}} 
\@tabular 
} 
% TABBING 
% 
\tabbingsep \labelsep 
 
% MINIPAGE 
% 
\skip\@mpfootins = \skip\footins 
 
% FRAMEBOX 
% 
% Space left between box and text by \fbox and \framebox. 
\setlength\fboxsep{3pt} 
% Width of rules in box made by \fbox and \framebox. 
\setlength\fboxrule{\arrayrulewidth}
 
%  **************************************** 
%  *         TITLE AND ABSTRACT           * 
%  **************************************** 
% 
\def\maketitle{\par 
 \begingroup 
%  \def\thefootnote{\fnsymbol{footnote}} 
%  \def\@makefnmark{\hbox to 0pt{$^{\@thefnmark}$\hss}} 
  \if@twocolumn 
   \twocolumn[\vspace*{17pt}\@maketitle]  
  \else 
   \newpage 
   \global\@topnum\z@ 
   \@maketitle 
  \fi 
  \thispagestyle{titlepage} 
%  \@thanks 
 \endgroup 
% \setcounter{footnote}{0} 
 \let\maketitle\relax 
 \let\@maketitle\relax 
 \gdef\@author{} 
 \gdef\@title{} 
% \gdef\@thanks{} 
 \let\thanks\relax 
} 
\def\and{\end{author@tabular}\vskip 3pt\par 
 \begin{author@tabular}[t]{@{}l@{}}} 
\def\@maketitle{\newpage 
 \vspace*{7pt} 
 {\raggedright \sloppy 
  {\reset@font\huge \bf \@title \par} 
  \vskip 23pt 
  {\reset@font\LARGE 
   \begin{author@tabular}[t]{@{}l@{}}\@author 
   \end{author@tabular}\par} 
  \vskip 22pt 
 } 
 \par\noindent 
 \ifGJ@mreferee
 \else
  {\reset@font\small \@date \par} 
   \vskip 22pt
 \fi 
} 
%  abstract/summary
\def\summary{\if@twocolumn 
   \start@GJbox\@summary 
 \else 
   \@summary 
 \fi} 
\def\endsummary{\if@twocolumn 
   \endlist\finish@GJbox 
 \else 
  \endlist 
 \fi} 
\def\abstract{\if@twocolumn 
   \start@GJbox\@summary 
 \else 
   \@summary 
 \fi} 
\def\endabstract{\if@twocolumn 
   \endlist\finish@GJbox 
 \else 
  \endlist 
 \fi} 
\def\left@zm{10.5pc}
\ifGJ@mreferee 
  \def\left@zm{1.5pc}
\fi
\ifGJ@referee 
  \def\left@zm{1.5pc} 
\fi
\def\@summary{\list{}{% 
    \labelwidth\z@ \labelsep\z@  
    \leftmargin\left@zm \rightmargin\z@ 
    \parsep 0pt plus 1pt}\item[]\reset@font
    \normalsize{\bf SUMMARY}\\\reset@font\large} % BK GJI 1.1 

% 
\newif\ifGJ@keywords 
\def\keywords{\if@twocolumn 
  \start@GJbox\@keywords 
 \else 
  \@keywords 
 \fi 
} 
\def\@keywords{\list{}{% 
    \labelwidth\z@ \labelsep\z@  
    \leftmargin\left@zm \rightmargin\z@ 
    \parsep 0pt plus 1pt}\item[]\reset@font\large{\bf Key words: }} 
\def\endkeywords{\if@twocolumn 
  \endlist\addvspace{37pt plus 0.5\baselineskip}\finish@GJbox  
 \else 
  \endlist 
 \fi 
 \@thanks 
 \gdef\@thanks{} 
 \GJ@keywordstrue 
} 
\def\nokeywords{\ifGJ@keywords\else 
 \if@twocolumn \start@GJbox\addvspace{37pt}\finish@GJbox \fi 
 \@thanks 
 \gdef\@thanks{}\fi 
} 
 
\def\author@tabular{\def\@halignto{}\@authortable} 
\let\endauthor@tabular=\endtabular 
% 
\def\author@tabcrone{\vspace{2pt}{\ifnum0=`}\fi% 
  \@xtabularcr[-7pt]\reset@font\small\it 
  \let\\=\author@tabcrtwo\ignorespaces} 
% 
\def\author@tabcrtwo{{\ifnum0=`}\fi\@xtabularcr[-7pt]\reset@font\small\it 
 \let\\=\author@tabcrtwo\ignorespaces} 
% 
\def\@authortable{\leavevmode \hbox \bgroup $\let\@acol\@tabacol 
 \let\@classz\@tabclassz \let\@classiv\@tabclassiv 
 \let\\=\author@tabcrone\ignorespaces \@tabarray} 
 
\def\start@GJbox{\@next\@currbox\@freelist{}{}% 
 \global\setbox\@currbox 
 \vbox\bgroup 
  \hsize \textwidth 
  \@parboxrestore 
} 
\def\finish@GJbox{\par\vskip -\dbltextfloatsep 
  \egroup 
  \global\count\@currbox\tw@ 
  \global\@dbltopnum\@ne 
  \global\@dbltoproom\maxdimen\@addtodblcol 
  \global\vsize\@colht 
  \global\@colroom\@colht 
} 
 
%  **************************************** 
%  *            PAGE STYLES               * 
%  **************************************** 
% 
\mark{{}{}} 
\gdef\@author{\mbox{}} 
\def\author{\@ifnextchar [{\@authortwo}{\@authorone}} 
\def\@authortwo[#1]#2{\gdef\@author{#2}\gdef\@shortauthor{#1}} 
\def\@authorone#1{\gdef\@author{#1}\gdef\@shortauthor{#1}} 
%\def\shortauthor#1{\gdef\@shortauthor{#1}} 
\gdef\@shortauthor{} 
\gdef\@title{\mbox{}} 
\def\title{\@ifnextchar [{\@titletwo}{\@titleone}} 
\def\@titletwo[#1]#2{\gdef\@title{#2}\gdef\@shorttitle{#1}} 
\def\@titleone#1{\gdef\@title{#1}\gdef\@shorttitle{#1}} 
%\def\shorttitle#1{\gdef\@shorttitle{#1}} 
\gdef\@shorttitle{} 
\def\volume#1{\gdef\@volume{#1}} 
\gdef\@volume{000} 
\def\tod@y{\today} 
\def\microfiche#1{\gdef\@microfiche{#1}} 
\gdef\@microfiche{} 
\def\pubyear#1{\gdef\@pubyear{#1}} 
\gdef\@pubyear{0000} 
\def\pagerange#1{\gdef\@pagerange{#1}} 
\gdef\@pagerange{000--000} 
\def\journal#1{\gdef\@journal{#1}} 
\def\plate#1{\gdef\pl@te{#1}} 
\ifx\CUP@mtlplain@loaded\undefined 
  \def\@printed{\qquad Printed \tod@y} 
\else 
  \def\@printed{} 
\fi 

\gdef\@journal{ } 
\ifGJ@referee
 \def\@journal{submitted to {\it Geophys. J. Int.} }
\fi
\ifGJ@mreferee
 \def\@journal{submitted to {\it Geophys. J. Int.}}
\fi
% 
% \ps@empty and \ps@plain defined in LATEX.TEX 
% 
\def\@ddell{} 
\if@fasttrack 
  \def\@ddell{F} 
\fi 
\def\ps@headings{\let\@mkboth\markboth 
\ifGJ@mreferee
 \def\@oddhead{\normalsize\hfill{\it\@shorttitle}\hspace{1.5em}\rm\@ddell\thepage} 
 \def\@oddfoot{} 
 \def\@evenhead{\normalsize\@ddell\thepage\hspace{1.5em}\it\@shortauthor\hfill} 
 \def\@evenfoot{} 
\else
 \def\@oddhead{\Large\hfill{\it\@shorttitle}\hspace{1.5em}\rm\@ddell\thepage} 
 \def\@oddfoot{} 
 \def\@evenhead{\Large\@ddell\thepage\hspace{1.5em}\it\@shortauthor\hfill} 
 \def\@evenfoot{} 
\fi
 \def\sectionmark##1{\markboth{##1}{}} 
 \def\subsectionmark##1{\markright{##1}}} 
\def\ps@myheadings{\let\@mkboth\@gobbletwo 
\ifGJ@mreferee
 \def\@oddhead{\normalsize\hfill\it\rightmark\hspace{1.5em}\rm\@ddell\thepage} 
 \def\@oddfoot{} 
 \def\@evenhead{\normalsize\@ddell\thepage\hspace{1.5em}\it\leftmark\hfill} 
 \def\@evenfoot{} 
\else
 \def\@oddhead{\Large\hfill\it\rightmark\hspace{1.5em}\rm\@ddell\thepage} 
 \def\@oddfoot{} 
 \def\@evenhead{\Large\@ddell\thepage\hspace{1.5em}\it\leftmark\hfill} 
 \def\@evenfoot{} 
\fi
 \def\sectionmark##1{} 
 \def\subsectionmark##1{}} 
\def\ps@titlepage{\let\@mkboth\@gobbletwo 
 \def\@oddhead{\footnotesize\@journal\hfill} 
%  \raisebox{0pt}[8pt][1pt]{\parbox{0.7\textwidth}{\@journal\\ 
%  \copyright\ \number\year\ BCS}} \hfill} 
 \def\@oddfoot{} 
 \def\@evenhead{\footnotesize\@journal\hfill} 
%  \raisebox{0pt}[8pt][1pt]{\parbox{0.7\textwidth}{\@journal\\ 
%  \copyright\ \number\year\ BCS}} \hfill} 
 \def\@evenfoot{} 
 \def\sectionmark##1{} 
 \def\subsectionmark##1{}} 
\def\ps@plate{\let\@mkboth\@gobbletwo 
  \def\@oddhead{\footnotesize\pl@te\hfill}\def\@oddfoot{} 
  \def\@evenhead{\footnotesize\pl@te\hfill}\def\@evenfoot{} 
} 
 
%  **************************************** 
%  *         TABLE OF CONTENTS, ETC.      * 
%  **************************************** 
% 
\def\@pnumwidth{1.55em} 
\def\@tocrmarg {2.55em} 
\def\@dotsep{4.5} 
\setcounter{tocdepth}{2} 
\def\@undottedtocline#1#2#3#4#5{\ifnum #1>\c@tocdepth 
 \else 
  \vskip \z@ plus .2pt 
  {\hangindent #2\relax 
   \rightskip \@tocrmarg \parfillskip -\rightskip 
   \parindent #2\relax \@afterindenttrue 
   \interlinepenalty\@M \leavevmode 
   \@tempdima #3\relax #4\nobreak \hfill \nobreak 
   \hbox to\@pnumwidth{\hfil\reset@font\rm \ }\par}\fi} 
\def\tableofcontents{\@restonecolfalse 
 \if@twocolumn\@restonecoltrue\onecolumn\fi 
 \section*{CONTENTS} \@starttoc{toc} 
 \if@restonecol\twocolumn\fi \par\vspace{12pt}} 
\def\l@part#1#2{\addpenalty{-\@highpenalty} 
 \addvspace{2.25em plus 1pt} 
 \begingroup 
  \parindent \z@ \rightskip \@pnumwidth 
  \parfillskip -\@pnumwidth 
  {\reset@font\normalsize\rm 
   \leavevmode \hspace*{3pc} 
   #1\hfil \hbox to\@pnumwidth{\hss \ }}\par 
   \nobreak \global\@nobreaktrue 
   \everypar{\global\@nobreakfalse\everypar{}}\endgroup} 
\def\l@section#1#2{\addpenalty{\@secpenalty} 
 \@tempdima 1.5em 
 \begingroup 
  \parindent \z@ \rightskip \@pnumwidth 
  \parfillskip -\@pnumwidth \reset@font\rm \leavevmode 
  \advance\leftskip\@tempdima \hskip -\leftskip 
  #1\nobreak\hfil \nobreak\hbox to\@pnumwidth{\hss \ }\par 
 \endgroup} 
\def\l@subsection{\@undottedtocline{2}{1.5em}{2.3em}} 
\def\l@subsubsection{\@undottedtocline{3}{3.8em}{3.2em}} 
\def\l@paragraph{\@undottedtocline{4}{7.0em}{4.1em}} 
\def\l@subparagraph{\@undottedtocline{5}{10em}{5em}} 
\def\listoffigures{\@restonecolfalse 
 \if@twocolumn\@restonecoltrue\onecolumn\fi 
 \section*{LIST OF FIGURES\@mkboth{LIST OF FIGURES}{LIST OF FIGURES}} 
 \@starttoc{lof} \if@restonecol\twocolumn\fi} 
\def\l@figure{\@undottedtocline{1}{1.5em}{2.3em}} 
\def\listoftables{\@restonecolfalse 
 \if@twocolumn\@restonecoltrue\onecolumn\fi 
 \section*{LIST OF TABLES\@mkboth{LIST OF TABLES}{LIST OF TABLES}} 
 \@starttoc{lot} \if@restonecol\twocolumn\fi} 
\let\l@table\l@figure 
 
%  **************************************** 
%  *             BIBLIOGRAPHY             * 
%  **************************************** 
% 
 
\newcounter{dummy} 
\def\thebibliography#1{\section*{REFERENCES} 
 \addcontentsline{toc}{section}{REFERENCES} 
 \list{}{\usecounter{dummy}% 
         \labelwidth\z@ 
         \leftmargin 1.5em 
         \itemsep \z@ 
         \itemindent-\leftmargin} 
 \reset@font\small 
 \parindent\z@ 
 \parskip\z@ plus .1pt\relax 
 \def\newblock{\hskip .11em plus .33em minus .07em} 
 \sloppy\clubpenalty4000\widowpenalty4000 
 \sfcode`\.=1000\relax 
} 
\let\endthebibliography=\endlist 
\def\@biblabel#1{\hspace*{\labelsep}[#1]} 

% For AGU compatability use a {references} environment
% It is not necessary to use the thebibliography environment.  This style
% defines a references environment that simply sets off the list of references
% and adjusts spacing parameters.  It is not possible to use \bibitems
% within the references environment.  (The \reference command cannot be
% shortened to \ref; \ref is used in vanilla LaTeX for cross-referencing.)
%
%    \begin{references}
%    \reference <bibliographic data>
%    \end{references}
%
%    \markcite{TEXT}
%
% Observe that the bibliographic data supplied by the author must also
% conform to the standards of the journal.  I have elected not to burden
% authors with tedious LaTeX commands which would require many curly braces
% to delimit the bibliographic fields because many of the journals I have
% targeted in this project have agreed to reduce typographic overhead (bolding,
% italicizing, etc.) in reference lists.  (See Abt's editorial in ApJ 357.)
% It is the responsibility of the author to get these fields in the proper
% order with the correct punctuation; the information will be typeset as is,
% i.e., in roman with no size or style changes.
%
% Macros are provided for many of the oft-referenced journals so that authors
% may use the LaTeX names rather than having to look up a particular journal's
% specific abbreviation.  In principle, all the journals should be using the
% same abbreviations, too, but I anticipate some changes in the specific
% abbreviations before a set is finally settled on.  As long as these macros
% are kept up to date, authors need not be concerned about such editorial
% changes.
%
% Corresponding with the non-BibTeX references environment and \reference
% command for the reference list is the \markcite command.  This can be
% used to identify in-text citations; in this case, the author provides
% the actual citation text as opposed to the symbolic tag.  The use of
% this is encouraged, even for journals where no special formatting of
% citation call-outs is required; it is crucial for the future use of
% on-line browsers.

\def\references{\section*{REFERENCES}
    \bgroup\parindent=\z@\parskip=\z@\reset@font\small
    \def\refpar{\par\hangindent=1em\hangafter=1}}
\def\endreferences{\refpar\egroup}
%%\let\endreferences=\endlist 

%% \small should probably be \scriptsize to match print version 


\def\reference{\relax\refpar}  % Generic: might be better...
\def\markcite#1{#1\relax}


 
%  **************************************** 
%  *              THE INDEX               * 
%  **************************************** 
% 
\newif\if@restonecol 
\def\theindex{\section*{INDEX} 
 \addcontentsline{toc}{section}{INDEX} 
 \reset@font\footnotesize \parindent\z@ \parskip\z@ plus .1pt\relax 
 \let\item\@idxitem} 
\def\@idxitem{\par\hangindent 1em} 
\def\subitem {\par\hangindent 1em \hspace*{1em}} 
\def\subsubitem{\par\hangindent 1em\hspace*{1em}} 
\def\endtheindex{\if@restonecol\onecolumn\else\clearpage\fi} 
\def\indexspace{\par\vskip\baselineskip \relax} 
 
%  **************************************** 
%  *             FOOTNOTES                * 
%  **************************************** 
% 
\def\footnoterule{\kern-3\p@ \hrule width 12pc height \z@ \kern 3\p@} 
 
\renewcommand{\thefootnote}{\mbox{{$\fnsymbol{footnote}$}}} 
 
\def\@fnsymbol#1{\ifcase#1\or \mbox{$\star$}\or \dagger\or \ddagger\or 
   \S \or \P \or \|\or **\or \dagger\dagger 
   \or \ddagger\ddagger\or \S\S\or \P\P\or \|\|\else *** 
   \fi\relax} 
 
\long\def\@makefntext#1{\parindent 1em\noindent 
  $^{\@thefnmark}$\hspace{4pt}#1} 
 
%\long\def\@makefntext#1{\@setpar{\@@par\@tempdima \hsize 
% \advance\@tempdima-18pt\parshape \@ne 18pt\@tempdima}\par 
% \noindent \hbox to \z@{\hss$^{\@thefnmark}$\hspace{0.5em}}#1} 
 
%  **************************************** 
%  *         FIGURES AND TABLES           * 
%  **************************************** 
% 
\setcounter{topnumber}{2} 
\renewcommand\topfraction{.9} 
\setcounter{bottomnumber}{1} 
\renewcommand\bottomfraction{.9} 
\setcounter{totalnumber}{3} 
\renewcommand\textfraction{.1} 
\renewcommand\floatpagefraction{.9} 
\setcounter{dbltopnumber}{2} 
\renewcommand\dbltopfraction{.9} 
\renewcommand\dblfloatpagefraction{.9} 
 
\newif\ifst@rredfloat 
\newcounter{table} 
\def\thetable{\@arabic\c@table} 
\def\fps@table{tbp} 
\def\ftype@table{1} 
\def\fstyle@table{\reset@font\small\rm} 
\def\ext@table{lot} 
\def\fnum@table{Table \thetable} 
\def\table{\let\@makecaption=\GJ@maketablecaption% 
  \@float{table}\fstyle@table} 
\let\endtable\end@float 
\@namedef{table*}{\st@rredfloattrue\let\@makecaption=\GJ@maketablecaption% 
  \@dblfloat{table}\centering\fstyle@table} 
\@namedef{endtable*}{\end@dblfloat\st@rredfloatfalse} 
 
\newcounter{figure} 
\def\thefigure{\@arabic\c@figure} 
\def\fps@figure{tbp} 
\def\ftype@figure{2} 
\def\fstyle@figure{\reset@font\small\rm} 
\def\ext@figure{lof} 
\def\fnum@figure{Figure \thefigure} 
\def\figure{\let\@makecaption=\GJ@makefigurecaption% 
  \@float{figure}\fstyle@figure} 
\let\endfigure\end@float 
\@namedef{figure*}{\st@rredfloattrue\let\@makecaption=\GJ@makefigurecaption% 
  \@dblfloat{figure}\centering\fstyle@figure} 
\@namedef{endfigure*}{\end@dblfloat\st@rredfloatfalse} 
 
\long\def\GJ@makefigurecaption#1#2{\vskip 6pt 
 \setbox\@tempboxa\hbox{\reset@font\small{\bf #1.} #2} 
 \ifdim \wd\@tempboxa >\hsize 
  \reset@font\small{\bf #1.} #2\par 
 \else 
   \ifst@rredfloat 
     \hbox to\hsize{\hfill\box\@tempboxa\hfill} 
   \else 
     \hbox to\hsize{\box\@tempboxa\hfill} 
   \fi 
 \fi 
 \vskip 6pt 
} 
 
\long\def\GJ@maketablecaption#1#2{\vskip 6pt 
 \setbox\@tempboxa\hbox{\reset@font\small{\bf #1.} #2} 
 \ifdim \wd\@tempboxa >\hsize 
   \reset@font\small{\bf #1.} #2\par 
 \else 
   \hbox to\hsize{\box\@tempboxa\hfill} 
 \fi 
 \vskip 6pt 
} 
 
%\def\caption{\refstepcounter\@captype \@dblarg{\@caption\@captype}} 
\def\caption{\@ifstar{\GJ@caption\@captype}% 
 {\refstepcounter\@captype \@dblarg{\@caption\@captype}}% 
} 
\long\def\GJ@caption#1#2{%\par 
 \begingroup 
  \@parboxrestore 
  \reset@font\normalsize 
  \@makecaption{\csname fnum@#1\endcsname}{\ignorespaces #2}\par 
 \endgroup} 
 
%  **************************************** 
%  *           Specials               * 
%  **************************************** 
% 
% 
%  The following code is taken from authordate.sty, by  David Rhead, of 
%  the Cripps Computing Centre (February 1990), with some modifications 
%  by ggw.  It is intended for documents 
%  that use the author-date citation system.  When used in conjunction with 
%  authordate1.bst, ... , authordate4.bst, it allows citations of the 
%  form \shortcite{bloggs-60} as well as the usual
%  \cite{bloggs-60}. 
%  When 
%       ... \cite{bloggs-60} ... \shortcite{bloggs-60} ... 
%  appears in the input file, 
%       ... (Bloggs, 1960) ... (1960) ... 
%  appears in the final document. 
% 
%  Code taken from apalike.sty. 
% 
\ifGJ@extra
\else
\def\@cite#1#2{(#1\if@tempswa , #2\fi)} 
\def\@cited#1#2{#1\if@tempswa (#2)\fi}
\def\@biblabel#1{} 
 
% Set length of hanging indentation for bibliography entries 
% 
\newlength{\bibhang} 
\setlength{\bibhang}{1em} 
\fi 
%  Define default title and page-head text. 
% 
\def\bibtitle{References} 
\def\bibheadtitle{REFERENCES} 
 
%  Code taken from aaai.sty. 
% 
% don't box citations, add space between multiple citations, separate with ; 
\def\@citex[#1]#2{\if@filesw\immediate\write\@auxout{\string\citation{#2}}\fi 
  \def\@citea{}\@cite{\@for\@citeb:=#2\do 
    {\@citea\def\@citea{; }\@ifundefined 
       {b@\@citeb}{@warning 
       {Citation `\@citeb' on page \thepage \space undefined}}% 
{\csname b@\@citeb\endcsname}}}{#1}} 
% Allow short (name-less) citations, when used in 
% conjunction with a bibliography style that creates labels like 
%       \citename{<names> }<year> 
\let\@internalcite\cite
\def\cite{\def\citename##1{##1}\@internalcite} 
\def\shortcite{\def\citename##1{}\@internalcite}

 
\def\bsp{\medskip\reset@font\small\noindent This paper has been produced  
 using the Blackwell Scientific Publications GJI \LaTeX2e\ class file.} 
 
 
%  **************************************** 
%  *                FLEQN                 * 
%  **************************************** 
% 
%% 
% FLEQN DOCUMENT-STYLE OPTION - released 9 March 1987 
%    for LaTeX version 2.09 
% Copyright (C) 1987 by Leslie Lamport 
% 
% The following commands make the displayed math environments flush left, 
% with an indentation of \mathindent from the prevailing left margin. 
 
\def\[{\relax\ifmmode\@badmath\else\begin{trivlist}\item[]\leavevmode 
  \hbox to\linewidth\bgroup$ 
  \displaystyle 
  \hskip\mathindent\bgroup\fi} 
 
\def\]{\relax\ifmmode \egroup $\hfil 
       \egroup \end{trivlist}\else \@badmath \fi} 
 
\def\equation{\refstepcounter{equation}\trivlist \item[]\leavevmode 
  \hbox to\linewidth\bgroup $ 
  \displaystyle 
\hskip\mathindent} 
 
 
\def\endequation{$\hfil 
           \displaywidth\linewidth\@eqnnum\egroup \endtrivlist} 
 
\def\eqnarray{\stepcounter{equation}\let\@currentlabel=\theequation 
\global\@eqnswtrue 
\global\@eqcnt\z@\tabskip\mathindent\let\\=\@eqncr 
\abovedisplayskip\topsep\ifvmode\advance\abovedisplayskip\partopsep\fi 
\belowdisplayskip\abovedisplayskip 
\belowdisplayshortskip\abovedisplayskip 
\abovedisplayshortskip\abovedisplayskip 
$$\halign 
to \linewidth\bgroup\@eqnsel\hskip\@centering$\displaystyle\tabskip\z@ 
  {##}$&\global\@eqcnt\@ne \hskip 2\arraycolsep \hfil${##}$\hfil 
  &\global\@eqcnt\tw@ \hskip 2\arraycolsep $\displaystyle{##}$\hfil 
   \tabskip\@centering&\llap{##}\tabskip\z@\cr} 
 
\def\endeqnarray{\@@eqncr\egroup 
 \global\advance\c@equation\m@ne$$\global\@ignoretrue} 
 
\newdimen\mathindent 
\mathindent = \z@ 
 
\ifGJ@mreferee 
 \setlength {\mathindent}  {2 mm}
\fi
 
%  **************************************** 
%  *            MISCELLANEOUS             * 
%  **************************************** 
% 
% DATE 
\def\today{\number\day\ \ifcase\month\or 
  January\or February\or March\or April\or May\or June\or 
  July\or August\or September\or October\or November\or December 
 \fi \ \number\year} 
 
%  **************************************** 
%  *           INITIALIZATION             * 
%  **************************************** 
% 
\flushbottom 
\ps@headings 
\ifGJ@galley 
 \raggedbottom 
 \ps@empty 
\fi 
\ifGJ@referee 
 \raggedbottom 
\fi 
\pagenumbering{arabic} 
\if@twocolumn 
 \twocolumn 
 \sloppy 
\else 
 \onecolumn 
\fi 
\frenchspacing 
 

%%%  ACKNOWLEDGEMENTS

\newenvironment{acknowledgments}{\section*{ACKNOWLEDGMENTS}}{}%
\newenvironment{acknowledgment}{\section*{ACKNOWLEDGMENT}}{}%
   

%%%  COMMON JOURNALS

\newcommand{\areps}{\textit{Ann.\ Rev.\ Earth Planet\ Sci.}, }
\newcommand{\bssa}{\textit{Bull.\ seism.\ Soc.\ Am.}, }
\newcommand{\eos}{\textit{Eos Trans.\ AGU}, } 
\newcommand{\eps}{\textit{Earth~Planets~Space},}
\newcommand{\epsl}{\textit{Earth~Planet.\ Sci.\ Lett.}, }
\newcommand{\gca}{\textit{Geochim.\ Cosmochim.\ Acta}, }
\newcommand{\geo}{\textit{Geology}, }
\newcommand{\geop}{\textit{Geophysics}, }
\newcommand{\gji}{\textit{Geophys.\ J\ Int.}, }
\newcommand{\gjras}{\textit{Geophys.\ J.\ R.\ astr.\ Soc.}, }
\newcommand{\grl}{\textit{Geophys.\ Res.\ Lett.}, }
\newcommand{\gsab}{\textit{Geol.\ Soc.\ Am.\ Bull.},}
\newcommand{\gs}{\textit{Geophys.\ Surv.}, }
\newcommand{\jgr}{\textit{J.\ geophys.\ Res.}, }
\newcommand{\jpo}{\textit{J.\ Phys.\ Oceanogr.}, } 
\newcommand{\jseis}{\textit{J.~Seismology}, }
\newcommand{\mnras}{\textit{Mon.\ Not.\ R.\ astr.\ Soc.}, } 
\newcommand{\pag}{\textit{Pure Appl.\ Geophys.},}
\newcommand{\pepi}{\textit{Phys.\ Earth Planet.\ Inter.}, }
\newcommand{\rg}{\textit{Rev.\ Geophys.}, }
\newcommand{\tecto}{\textit{Tectonophysics}, }

\ifGJ@extra
\def\refname{REFERENCES}
\usepackage{gji_extra}
\fi
% end of GJI.cls 
