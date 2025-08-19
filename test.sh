# python eval.py --profile profile/eval/s_i1_base.txt --mode face --weighted --csv resnet_results.csv
# python eval.py --profile profile/eval/s_i1_new.txt --mode face --weighted --csv mamba_results.csv
# python compare.py --new mamba_results.csv --base resnet_results.csv --labels "Mamba,ResNet"

python eval.py --profile profile/eval/s_i2_base.txt --mode face --weighted --csv resnet_results2.csv
python eval.py --profile profile/eval/s_i2_new.txt --mode face --weighted --csv mamba_results2.csv
python compare.py --new mamba_results2.csv --base resnet_results2.csv --labels "Mamba,ResNet"