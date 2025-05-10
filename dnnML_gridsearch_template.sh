#!/bin/bash
source /usr/local/root/6.32.04/bin/thisroot.sh
cd /home/xiaoxue/code/CMSTool/QCDMachine_Leaning
InputFolder="/storage/xiaoxue/code/data/Features_250510/"
InputFile="QCD_Pt-15to7000_Flat_herwig7-Run3Summer22EE_SpinOFF_CMSPartonReco_jinitpt_700_7000_j2pt_55_7000"
outputFolder="GridSearch_Fourlabel_reco"
rm -r ${outputFolder}
mkdir -p ${outputFolder}
cd ${outputFolder}
mkdir -p "bestmodel"
mkdir -p "confusion_matrix"
mkdir -p "loss_epoch"
mkdir -p "models"
mkdir -p "roc"
cd ..

python3 dnn_load_train_ak8_gridserach.py --sample_path0 "${InputFolder}${InputFile}/*.root" --sample_path1 ${outputFolder} --entries -1 --mode Train