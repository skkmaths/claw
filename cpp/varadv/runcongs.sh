NCELL=$1
FILE=error.txt
rm -rf $FILE && touch $FILE
M_PI=$(echo "scale=20; 4*a(1)" | bc -l)
for ncell in $NCELL
do 
   echo "ncell = $ncell"
   ./run -nx $ncell -ny $ncell -Tf $(echo "2.0 * $M_PI" | bc -l) -cfl 0.4 -save_freq 0 -scheme so >log.txt
   tail -n 1 log.txt
   tail -n 1 log.txt >> $FILE
done
echo "Wrote file $FILE"