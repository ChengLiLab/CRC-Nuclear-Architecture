for i in `ls ../HIC*`;do `nohup python /lustre/user/liclab/ganjb/tools/hicpro/HiC-Pro_2.11.1/bin/utils/sparseToDense.py -b 100k_abs.bed $i --perchr > ${i##*/}.log 2>&1 &`;done
