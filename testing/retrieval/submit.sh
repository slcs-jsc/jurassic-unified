if [ "$#" -ne 1 ];  then
  echo "Illegal number of parameters"
  exit 1
fi

if [ "$1" = "unified" ] || [ "$1" = "reference" ]; then
    echo "starting script for parameter: " $1
else
    echo "wrong input parameter"
    exit 1
fi

if [ "$1" = "unified" ]; then
  src="../../jurassic-gpu/src"
  name=$src/sca_retrieval
fi

if [ "$1" = "reference" ]; then
  src="../../../../../experiments/ref_repos/jurassic-scatter/src"
  name=$src/retrieval
fi

make -C $src

echo "rm out-$1"
rm out-$1
echo "rm err-$1"
rm err-$1

read -n 1 -s -r -p "Press any key to continue"
echo

sbatch jurun-cloud.sh $name

watch -n 1 squeue -u pozgaj1

echo "mv out out-"${1}
mv out out-$1
echo "mv err err-"${1}
mv err err-$1

watch -n 1 squeue -u pozgaj1
