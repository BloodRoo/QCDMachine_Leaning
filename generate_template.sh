#!/bin/bash
source /usr/local/root/6.32.04/bin/thisroot.sh
BaseFolder="/storage/xiaoxue/code/data/"
cd /home/xiaoxue/code/CMSTool/QCDMachine_Leaning
#################### 输入 ####################
Version="250510"
InputFolder="${BaseFolder}herwig/"
# InputFolder="${BaseFolder}CMSMC/AK8/AK8_MC/"
jinit_ptlow="700"
jinit_pthigh="7000" 
j2_ptlow="130"
j2_pthigh="7000"
SampleType="PrivateMC"
# SampleType="PrivateMCParton"
# SampleType="CMSPartonReco"
# SampleType="CMSData"
# SampleType="CMSMC"
InputFiles=(
              "QCD_Pt-15to7000_Flat_herwig7-Run3Summer22"
              "QCD_Pt-15to7000_Flat_herwig7-Run3Summer22EE"
              "QCD_Pt-15to7000_Flat_herwig7-Run3Summer22_SpinOFF"
)
#################### 输入 ####################

compile generate_sample_${Version}.cpp
for InputFile in "${InputFiles[@]}"; do
    target_dir="${BaseFolder}Machine_Learning/Features_${Version}/${InputFile}_${SampleType}_jinitpt_${jinit_ptlow}_${jinit_pthigh}_j2pt_${j2_ptlow}_${j2_pthigh}"
    [[ ! -d $target_dir ]] && mkdir -p $target_dir
    # ./generate_sample -i 0 -n 50 -k 50 -I ${InputFolder}${InputFile} -O /storage/xiaoxue/code/data/Machine_Learning/Features/${InputFile}_${SampleType}_jinitpt_${jinit_ptlow}_${jinit_pthigh}_j2pt_${j2_ptlow}_${j2_pthigh} --jinit_ptlow ${jinit_ptlow} --jinit_pthigh ${jinit_pthigh} --j2_ptlow ${j2_ptlow} --j2_pthigh ${j2_pthigh} --SampleType ${SampleType}
    run_parallel generate_sample -nchunks 50 -nparts 50 -opt "-I ${InputFolder}${InputFile} -O ${target_dir} --jinit_ptlow ${jinit_ptlow} --jinit_pthigh ${jinit_pthigh} --j2_ptlow ${j2_ptlow} --j2_pthigh ${j2_pthigh} --SampleType ${SampleType}"
done