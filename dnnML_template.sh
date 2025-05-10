#!/bin/bash
source /usr/local/root/6.32.04/bin/thisroot.sh
InputFolder="/storage/xiaoxue/code/data/Features_250510/"
Model_Name="/storage/shuangyuan/code/analysis_spin/Machine_learning/GridSearch_Fourlabel_new/models/Foulabel_hidden256_128_32_lr0.0005_l20_batch256_dropout0.0"
# Model_Name="/home/xiaoxue/code/spin/ways/Machine_Learning/GridSearch_Fourlabel_reco/models/Foulabel_hidden_256_128_128_lr0.0005_l20_batch256_dropout0.0"
InputFiles=(
            "QCD_PT-470to600_pythia8-Run3Summer22_CMSMC_jinitpt_700_7000_j2pt_55_7000"
            "QCD_PT-600to800_pythia8-Run3Summer22_CMSMC_jinitpt_700_7000_j2pt_55_7000"
            "QCD_PT-800to1000_pythia8-Run3Summer22_CMSMC_jinitpt_700_7000_j2pt_55_7000"
            "QCD_PT-1000to1400_pythia8-Run3Summer22_CMSMC_jinitpt_700_7000_j2pt_55_7000"
            "QCD_PT-1400to1800_pythia8-Run3Summer22_CMSMC_jinitpt_700_7000_j2pt_55_7000"
            "QCD_PT-1800to2400_pythia8-Run3Summer22_CMSMC_jinitpt_700_7000_j2pt_55_7000"
            "QCD_PT-2400to3200_pythia8-Run3Summer22_CMSMC_jinitpt_700_7000_j2pt_55_7000"
            "QCD_PT-2400to3200_pythia8-Run3Summer22_CMSMC_jinitpt_700_7000_j2pt_55_7000"
            "QCD_PT-3200toInf_pythia8-Run3Summer22_CMSMC_jinitpt_700_7000_j2pt_55_7000"
)

for InputFile in "${InputFiles[@]}"; do
    outputFolder=${InputFolder}${InputFile}
    python3 dnn_load_train_ak8.py --sample_path0 "${outputFolder}/*.root" --entries -1 --mode Prediction --model_path ${Model_Name}
done