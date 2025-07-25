# ./LSD-Gdata_mt profile/s_i1/profile1.txt
# ./LSD-Gdata_mt profile/s_i1/profile2.txt

# python train.py
# cp -r out/* model/s_i1

# ./LSD-denoising_mt profile/s_i1/s_i1.txt

./LSD-Gdata_mt profile/s_i2/profile1.txt
./LSD-Gdata_mt profile/s_i2/profile2.txt
python train.py
cp -r out model/s_i2