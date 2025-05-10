import sys, os
path = "X:/storage/ylye/QCDAnalysis/OfflineExamples/python" 
path = "python" 
sys.path.append(path)
from keras.optimizers import Adam
import tensorflow as tf
import pandas as pd
import uproot
import numpy as np
import argparse
from sklearn.model_selection import train_test_split
from tqdm import tqdm
import os
import math
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler, MinMaxScaler
import time
import dnn
from tensorflow.keras.models import load_model # type: ignore

def main(args):
    sample_path0 = args.sample_path0
    sample_path1 = args.sample_path1
    model_path = args.model_path
    suffix =args.suffix
    mode= args.mode
    
    dir_path = os.path.dirname(sample_path0)
    basefolder_name = os.path.basename(dir_path)+"_FourLabel1/"    
    select_item = ["nparticles2","ntracks2","pTD2","sigma12","sigma22","sigma2","nparticles3","ntracks3","pTD3","sigma13","sigma23","sigma3","nparticles4","ntracks4","pTD4","sigma14","sigma24","sigma4","z1","z2","deltaR2","kt2"]
    
    entrynum_iter=100000
    #entrynum_iter=10
    tree,nentries =  dnn.LoadTree(sample_paths=[sample_path0])
    if nentries>args.entries and args.entries!=-1:
        nentries = args.entries
    iter_num = math.ceil(nentries / entrynum_iter)
    model = tf.keras.Sequential()
    X_data1_train_full=np.empty((0, 0))
    X_data1_val_full=np.empty((0, 0))
    Y_data1_train_full=np.empty((0, 0))
    Y_data1_val_full=np.empty((0, 0))
    train_val_index_full=np.empty((0, 0))
    for iter in range(iter_num):
        X_data, branch_to_index, index_to_branch,branch_names,branch_types = dnn.LoadROOTFile(sample_paths=[sample_path0],entries=entrynum_iter,entry_begin=iter*entrynum_iter)
        if mode == "Train":
            folder_name = '/storage/shuangyuan/code/analysis_spin/Machine_learning/ML/Datasets/Train/'+basefolder_name
            if not os.path.exists(folder_name):
                os.makedirs(folder_name)
            match_index=["match2"]
            Y_data_index=["type2"]
            final_descri3_index=["final_descri3"]
            final_descri4_index=["final_descri4"]
            level_prefix=""
            level_suffix="_Hadron"
            #X_data_supp, *_ = dnn.LoadROOTFile(sample_paths=[sample_path1],entries=args.entries)
            X_data0,X_data_shape =  dnn.Extract_flat_features(X_data,branch_to_index,select_item,level_prefix,level_suffix,return_shape=True)
            # 创建一个StandardScaler对象
            scaler = StandardScaler()
            # 拟合并转换数据（计算均值并执行中心化）
            # X_data0 = scaler.fit_transform(X_data0)
            #X_data_supp0,X_data_supp_shape =  dnn.Extract_flat_features(X_data_supp,branch_to_index,select_item,level_prefix,level_suffix,return_shape=True)
            Y_data,*_ = dnn.Extract_flat_features(X_data,branch_to_index,Y_data_index,level_prefix,level_suffix,return_shape=True)
            match_data,*_ = dnn.Extract_flat_features(X_data,branch_to_index,match_index,prefix="",suffix="",return_shape=True)
            final_descri3_data,*_ = dnn.Extract_flat_features(X_data,branch_to_index,final_descri3_index,level_prefix,level_suffix,return_shape=True)
            final_descri4_data,*_ = dnn.Extract_flat_features(X_data,branch_to_index,final_descri4_index,level_prefix,level_suffix,return_shape=True)
            match_data = match_data[:,0]
            final_descri3_data = final_descri3_data[:,0]
            final_descri4_data = final_descri4_data[:,0]
            #Y_data0 = dnn.convert_to_one_hot_threeLabel(Y_data,match_data)
            Y_data0 = dnn.convert_to_one_hot_fourlabel3(Y_data,match_data,final_descri3_data,final_descri4_data)
            indices_type = dnn.balance_data_multiLabel(Y_data0)
            # indices_type = dnn.Selection(X_data0,indices_type)
            # indices_type = np.ones(len(indices_type), dtype=bool)
            X_data1 = X_data0[indices_type] # type: ignore
            Y_data1 = Y_data0[indices_type]
            select_index_train,select_index_val = dnn.custom_train_test_split_indices(len(X_data1),test_size=0.2,random_state=42)
            X_data1_train = X_data1[select_index_train]
            X_data1_val = X_data1[select_index_val]
            Y_data1_train = Y_data1[select_index_train]
            Y_data1_val = Y_data1[select_index_val]
            ###### X data full
            if iter == 0:
                X_data1_train_full=X_data1_train 
                X_data1_val_full=X_data1_val
                Y_data1_train_full=Y_data1_train
                Y_data1_val_full=Y_data1_val
            else:
                X_data1_train_full = np.concatenate((X_data1_train_full, X_data1_train), axis=0)
                X_data1_val_full = np.concatenate((X_data1_val_full, X_data1_val), axis=0)
                Y_data1_train_full = np.concatenate((Y_data1_train_full, Y_data1_train), axis=0)
                Y_data1_val_full = np.concatenate((Y_data1_val_full, Y_data1_val), axis=0)
            X_data1_train_full = scaler.fit_transform(X_data1_train_full)
            X_data1_val_full = scaler.transform(X_data1_val_full)
            ###### Add validation branch ##########
            train_val_index0 = dnn.modify_array_advanced(indices_type,select_index_train,select_index_val)
            train_val_index = dnn.ShapeAlign(train_val_index0,X_data_shape) 
            train_val_index_full = train_val_index_full + train_val_index if len(train_val_index_full) != 0 else train_val_index
            #train_val_index = dnn.ShapeAlign(train_val_index0,X_data_shape) +dnn.CreateAlignedShapeArr(X_data_supp_shape, -1)
            del X_data,X_data0, X_data1, Y_data, Y_data0, Y_data1, match_data, final_descri3_data, final_descri4_data, select_index_train, select_index_val, X_data1_train, X_data1_val, Y_data1_train, Y_data1_val, train_val_index0, train_val_index
        if mode == "Prediction":
            model_setting=os.path.basename(model_path)
            folder_name = '/storage/xiaoxue/code/data/Machine_Learning/rootfile_storage/score_'+model_setting+'/'+basefolder_name
            if not os.path.exists(folder_name):
                os.makedirs(folder_name)
            model = load_model(model_path)
            Xdata_full = X_data
            level_prefixs = set()
            level_suffixs = set()
            for name in branch_names:
                if "Gen" in name and name !='GeneratorWeight':
                    level_prefixs.add("Gen")
                if "Parton" in name and "_Parton" not in name:
                    level_prefixs.add("Parton")
                if "Reco" in name and name !='RecoBtag_value':
                    level_prefixs.add("Reco")
                if "_Hadron" in name:
                    level_suffixs.add("_Hadron")
            level_prefixs = list(level_prefixs)
            level_suffixs = list(level_suffixs)
            if len(level_prefixs) == 0:
                level_prefixs=[""]
            if len(level_suffixs) == 0:
                level_suffixs=[""]   
            for level_prefix in level_prefixs:
                for level_suffix in level_suffixs:
                    
                    X_data_predict,X_data_shape =  dnn.Extract_flat_features(Xdata_full,branch_to_index,select_item,level_prefix,level_suffix,return_shape=True)
                    full_predictions = model.predict(X_data_predict,batch_size=1024)
                    n_classes = full_predictions.shape[1]
                    full_predictions =  dnn.ShapeAlign(full_predictions.T,X_data_shape)
                    for i in range(n_classes):
                        branch_names.append(level_prefix+f'score{i}'+level_suffix) 
                        branch_types.append('vector<double>') 
                        Xdata_full = [row + [new_value] for row, new_value in zip(Xdata_full, full_predictions[i])]
            if len(args.train_val_index_full)>0:
                entry_end = (iter+1)*entrynum_iter if len(args.train_val_index_full) > (iter+1)*entrynum_iter else len(args.train_val_index_full)
                train_val_index = args.train_val_index_full[iter*entrynum_iter:entry_end]
                branch_names.append('validate'+level_suffix) 
                branch_types.append('vector<int>') 
                Xdata_full = [row + [new_value] for row, new_value in zip(Xdata_full, train_val_index)]
            print("Save to root file")
            ###### Save to root file ###############
            dnn.save_multiclass_predictions_to_root(Xdata_full,root_filename=folder_name+f"Chunk{iter}",branch_names=branch_names,branch_types=branch_types)
            del X_data, Xdata_full, X_data_predict, full_predictions
    if mode == "Train":
        ######################### model train : X_data3_train , X_data3_val , Y_data1_train , Y_data1_val #####################################
        hidden_units=[128, 64]
        # hidden_units=[256,128,64]
        learning_rate=0.0001
        l2_reg=0.0001
        file_suffix='herwig730_off_ak8_standardscaler_3'
        model = dnn.train_and_save_model_MultiLabel(X_train=X_data1_train_full,X_val=X_data1_val_full,Y_train=Y_data1_train_full,Y_val=Y_data1_val_full,hidden_units=hidden_units,learning_rate=learning_rate,l2_reg=l2_reg,model=model,file_suffix=file_suffix)
        #model.save("ML/model/Threelabel_"+suffix)
        model.save("ML/model/Fourlabel1_"+suffix)
        ######################## roc curve ##############################
        predictions_train = model.predict(X_data1_train_full,batch_size=256)
        predictions_val = model.predict(X_data1_val_full,batch_size=256)
        dnn.plot_roc_curve_threeLabel(Y_data1_train_full,predictions_train,folder_name+"train"+suffix)
        dnn.plot_roc_curve_threeLabel(Y_data1_val_full,predictions_val,folder_name+"val"+suffix)
        #args.model_path = "ML/model/Threelabel_"+suffix
        args.model_path = "ML/model/Fourlabel1_"+suffix
        args.mode = "Prediction"
        args.train_val_index_full = train_val_index_full
        del X_data1_train_full,X_data1_val_full,Y_data1_train_full,Y_data1_val_full
        main(args)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the binary classification model.")
    parser.add_argument('--sample_path0', type=str, default='', help='Path to the root file for training or prediction (typically spin-off samples).')
    parser.add_argument('--sample_path1', type=str, default='', help='Path to the root file for validation (typically spin-on samples).')

    parser.add_argument('--mode', type=str, default="Prediction", help='Mode of operation: "Train" for model training, "Prediction" for model prediction.')
    parser.add_argument('--model_path', type=str, default="", help='Path to the model for prediction (required if mode is "Prediction").')
    parser.add_argument('--entries', type=int, default=-1, help='Number of entries for training or prediction.')
    parser.add_argument('--suffix', type=str, default="", help='Suffix for the name of saved model and roc curve.')
    
    args = parser.parse_args()
    args.train_val_index_full =[]
    main(args)