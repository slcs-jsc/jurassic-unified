make clean -C ../../jurassic-gpu/src/
make -C ../../jurassic-gpu/src/

rm out err rad.tab

ls ../../jurassic-gpu/src/ -rtl | tail -2
# ls ../../src/ -rtl | tail -1

read -n 1 -s -r -p "Press any key to continue"
echo
sbatch jurun.sh ../../jurassic-gpu/src $m

watch -n 1 squeue -u pozgaj1

python3 diff.py generated_data/rad.org rad.tab
read -n 1 -s -r -p "Press any key to continue"
echo
