NCELL=$1
FILE=error.txt
rm -rf $FILE && touch $FILE
for ncell in $NCELL
do 
   echo "ncell = $ncell"
   python clw.py -ic smooth -Tf 0.13 -cfl 0.9 -nc $ncell -pde burger -numflux nt -compute_error yes -limit no \
          -plot_freq 0 -time_scheme euler>log.txt
   tail -n 1 log.txt
   tail -n 1 log.txt >> $FILE
done
echo "Wrote file $FILE"