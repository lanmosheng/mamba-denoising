# bash build.sh
# ./LSD-Gdata_mt profile/s_i1/profile1.txt
# ./LSD-Gdata_mt profile/s_i1/profile2.txt
# python train.py 40 40

python build_centers.py \
  --patch-root patches \
  --iou 0.75 --shuffle

python train.py
