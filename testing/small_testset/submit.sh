make -C ../../src/

rm out
n=$(< aux/submission_index)
m=$(( n + 1 ))
echo $m > aux/submission_index
echo last: $n, new: $m
ls ../../src/ -rtl | tail -3

read -n 1 -s -r -p "Press any key to continue"
echo
sbatch jurun-ice-785.sh ../../src $m

watch -n 1 squeue -u pozgaj1

python3 diff.py
read -n 1 -s -r -p "Press any key to continue"
echo
less out
