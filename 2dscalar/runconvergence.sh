NCELL=$1
FILE=error.txt
rm -rf $FILE && touch $FILE
for ncell in $NCELL
do 
   echo "ncell = $ncell"
   python clw2d.py -Tf 1.0 -ncellx $ncell -ncelly $ncell -compute_error yes \
          -plot_freq 0 -scheme mh -limit kappa -cfl 0.5 -save_freq 0 >log.txt
   tail -n 1 log.txt
   tail -n 1 log.txt >> $FILE
done
python plotrate.py
echo "Wrote file $FILE"