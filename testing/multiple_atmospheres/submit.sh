make clean -C ../../src/
make -C ../../src/

rm out err rad.tab
rm rad0.tab rad1.tab rad2.tab

ls ../../src/ -rtl | tail -3

read -n 1 -s -r -p "Press any key to continue"
echo
sbatch jurun.sh ../../src $m

watch -n 1 squeue -u pozgaj1

python3 diff.py generated_data/rad0.org rad0.tab
python3 diff.py generated_data/rad1.org rad1.tab
python3 diff.py generated_data/rad2.org rad2.tab
read -n 1 -s -r -p "Press any key to continue"
echo
