g++ -o greedymatching greedymatching.cpp
for i in new-*.txt;
do  echo "$i"; ./greedymatching    < "$i" > "$i-GreedyMatching"; done
