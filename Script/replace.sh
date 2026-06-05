for i in *.itp
do
sed -i 's/DD1/CC1/g' $i
sed -i 's/D1/C1/g' $i
sed -i 's/AA1/CC1/g' $i
sed -i 's/A1/C1/g' $i
done
