#! /bin/bash
set -e

mkdir -p ,,latex
for i in chapter0*.tm ; do
  TEXFILE=,,latex/${i%.tm}.tex
  TEMPTEXFILE=,,latex/${i%.tm}temp.tex
  if test ! -f $TEXFILE -o "(" $i -nt $TEXFILE ")"; then
    echo converting $i to $TEXFILE
    texmacs -c $i $TEMPTEXFILE --quit
    python strip_texmacs_boilerplate.py $TEMPTEXFILE $TEXFILE
  fi
  rm $TEMPTEXFILE
done

for i in *.fig; do
    fig2dev -Lpstex $i ,,latex/${i%.fig}.pstex
    fig2dev -Lpstex_t -p${i%.fig}.pstex $i ,,latex/${i%.fig}.pstex_t
done

cd ,,latex
export TEXINPUTS=$TEXINPUTS:..
cp ../multigrid.tex .
latex multigrid.tex || { less multigrid.log; exit 1; }
latex multigrid.tex
latex multigrid.tex
dvips multigrid.dvi -o multigrid.ps
ps2pdf multigrid.ps
