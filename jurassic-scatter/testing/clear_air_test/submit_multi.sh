if [ "$#" -ne 1 ];  then
  echo "Illegal number of parameters"
  exit 1
fi

if [ "$1" = "gpu" ] || [ "$1" = "plain" ] || [ "$1" = "scatter" ] || [ "$1" = "unified-sca" ] || [ "$1" = "unified-gpu" ]; then
    echo "starting script for parameter: " $1
else
    echo "wrong input parameter"
    exit 1
fi

echo $1 > last_submit_type

if [ "$1" = "unified-sca" ] || [ "$1" = "unified-gpu" ]; then
  src="../../jurassic-gpu/src"
fi

if [ "$1" = "gpu" ] || [ "$1" = "plain" ] || [ "$1" = "scatter" ]; then
  src="../../../reference_projects/jurassic-"$1"/src"
fi

make clean -C $src
make -C $src

echo "rm out-$1"
rm out-$1
rm err-$1

echo "rm rad-"${1}".tab"
rm rad-${1}.tab

read -n 1 -s -r -p "Press any key to continue"
echo

if [ "$1" = "unified-sca" ] || [ "$1" = "unified-gpu" ]; then
  sbatch jurun-unified.sh $1
else
  if [ "$1" = "gpu" ]; then
    sbatch jurun-reference-gpus.sh $src $1
  else
    sbatch jurun-reference-batch.sh $src $1
  fi
fi

watch -n 1 squeue -u pozgaj1

echo "mv out out-"${1}
mv out out-$1
echo "mv err err-"${1}
mv err err-$1
