NCELL=$1
FILE=error.txt
rm -rf $FILE && touch $FILE
M_PI=$(echo "scale=20; 4*a(1)" | bc -l)
printf "%-8s %-8s %-12s %-16s\n" \
"#cells" "dx" "l1error" "wallclock" > "$FILE"
for ncell in $NCELL
do 
   echo "ncell = $ncell"
   ./run -nx $ncell -ny $ncell -Tf 0.1 -cfl 0.5 -save_freq 0 -scheme fo -flux_type uw >log.txt
   tail -n 1 log.txt
   printf "%-8s %-8s %-12s %-16s\n" $(tail -n 1 log.txt) >> "$FILE"
done
echo "Wrote file $FILE"

# Note
#------------------------------------------
# to test EOC for linadv on [-1,1]X[-1,1]
# set Tf=1, advection_velocity to (1,1)
# keep cfl = 0.5 for obtain better EOC
# first try uw flux
#-----------------------------------------

# for Tf = 2*PI, use -Tf $(echo "2.0 * $M_PI" | bc -l) 