rm -v -r unified_library
mkdir -v unified_library
cp -v -r src/* unified_library

echo cd unified_library
cd unified_library
rm -v -r *.o
ls | grep -v "\." | xargs rm

cp -v ../src/Makefile .

declare -a structs=("aero_t" "atm_t" "ctl_t" "obs_t" "pos_t" "queue_item_t" "queue_t" "ret_t" "trans_table_t")

for val in ${structs[@]}; do
  echo substitute $val with jur_$val
  find ./ -type f -exec sed -i -e 's/'$val'/jur_'$val'/g' {} \;
done

make
