NCELL=$1
FILE=error.txt
rm -rf $FILE && touch $FILE
for ncell in $NCELL
do 
   echo "ncell = $ncell"
   python clw.py -ic smooth -Tf 2.0 -cfl 0.9 -nc $ncell -pde linear -numflux rusanov  -compute_error yes -theta 0.5 \
          -plot_freq 0 -time_scheme ssprk22>log.txt
   tail -n 1 log.txt
   tail -n 1 log.txt >> $FILE
done
echo "Wrote file $FILE"