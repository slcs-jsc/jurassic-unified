# src=../../../jurassic-gpu/src
src="../../../../../../experiments/ref_repos/jurassic-scatter/src"
make allclean -C $src
make -C $src

rm out
ls $src -rtl | tail -3

read -n 1 -s -r -p "Press any key to continue"
echo
sbatch jurun.sh $src

watch -n 1 squeue -u pozgaj1

python3 diff.py rad.org rad.tab
read -n 1 -s -r -p "Press any key to continue"
echo
less out
