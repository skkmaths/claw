NCELL=$1
FILE=error.txt
rm -rf $FILE && touch $FILE
make clean
make
for ncell in $NCELL
do 
   echo "ncell = $ncell"
   gmsh -2 strumesh.geo -setnumber lc $ncell -o mesh.msh
   ./run >log.txt
   tail -n 1 log.txt
   tail -n 1 log.txt >> $FILE
done
echo "Wrote file $FILE"