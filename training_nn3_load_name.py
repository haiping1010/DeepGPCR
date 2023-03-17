import numpy as np
import pandas as pd
import sys, os
from random import shuffle
import torch
import torch.nn as nn
from models.gat import GATNet
from models.gat_gcn import GAT_GCN
from models.gcn import GCNNet
from models.ginconv import GINConvNet
#from  read_smi_protein import *
from torch_geometric.data import InMemoryDataset, DataLoader
from torch_geometric import data as DATA
from utils import *

from sklearn.metrics import matthews_corrcoef


from sklearn import metrics

def acc(true, pred):

    return np.sum(true == pred) * 1.0 / len(true)

def aucJ(true_labels, predictions):

    fpr, tpr, thresholds = metrics.roc_curve(true_labels, predictions, pos_label=1)
    auc = metrics.auc(fpr,tpr)

    return auc




# training function at each epoch
def train(model, device, train_loader, optimizer, epoch, TRAIN_BATCH_SIZE):
    print('Training on {} samples...'.format(len(train_loader.dataset)))
    model.train()
    for batch_idx, data in enumerate(train_loader):
        data = data.to(device)
        optimizer.zero_grad()
        output = model(data,TRAIN_BATCH_SIZE,device)
        #output = model(data)
        loss = loss_fn(output, data.y.view(-1, 1).float().to(device))
        loss.backward()
        optimizer.step()
        if batch_idx % LOG_INTERVAL == 0:
            print('Train epoch: {} [{}/{} ({:.0f}%)]\tLoss: {:.6f}'.format(epoch,
                                                                           batch_idx * len(data.x),
                                                                           len(train_loader.dataset),
                                                                           100. * batch_idx / len(train_loader),
                                                                           loss.item()))

def predicting(model, device, loader,TRAIN_BATCH_SIZE):
    model.eval()
    total_preds = torch.Tensor()
    total_labels = torch.Tensor()
    names= []
    print('Make prediction for {} samples...'.format(len(loader.dataset)))
    with torch.no_grad():
        for data in loader:
            data = data.to(device)
            output = model(data,TRAIN_BATCH_SIZE,device)
            names=names+data.name
            total_preds = torch.cat((total_preds, output.cpu()), 0)
            total_labels = torch.cat((total_labels, data.y.view(-1, 1).cpu()), 0)
    return total_labels.numpy().flatten(),total_preds.numpy().flatten(),names


#modeling = [GINConvNet, GATNet, GAT_GCN, GCNNet][int(sys.argv[2])]
#model_st = modeling.__name__
modeling = GCNNet
model_st = modeling.__name__

import sys
sys.setrecursionlimit(100000)

#cuda_name = "cuda:0"
cuda_name = "cuda"
print (str(len(sys.argv))+'xxxx')
if len(sys.argv)>1:
    cuda_name = ["cuda:0","cuda:1"][int(sys.argv[1])]
print('cuda_name:', cuda_name)
#cuda_name = "cuda:0"
TRAIN_BATCH_SIZE = 50
LR = 0.0005
LOG_INTERVAL = 20
NUM_EPOCHS = 1

print('Learning rate: ', LR)
print('Epochs: ', NUM_EPOCHS)

import glob
from torch.utils.data.dataset import Dataset, ConcatDataset

dataset=TestbedDataset2(root='data1', dataset='L_P_test_neg')

from torch.utils.data import random_split
np.random.seed(0)
torch.manual_seed(0)

import torch
from torch.utils.data import Subset
from sklearn.model_selection import train_test_split

def train_val_dataset(dataset, val_split=0.25):
    train_idx, val_idx = train_test_split(list(range(len(dataset))), test_size=val_split)
    datasets = {}
    datasets['train'] = Subset(dataset, train_idx)
    datasets['val'] = Subset(dataset, val_idx)
    return datasets

dataset_pos=TestbedDataset2(root='data1', dataset='L_P_test_pos')

pos_num=dataset_pos.__len__()
neg_num=dataset.__len__()

print (dataset_pos.__len__(), dataset.__len__())
dataset_all=dataset_pos+dataset
print (torch.cuda.get_device_name(0))
print (torch.cuda.is_available())
# Main program: iterate over different datasets
for datasetxxxx in 'L':
        test_loader = DataLoader(dataset_all, batch_size=TRAIN_BATCH_SIZE, shuffle=True)
        print('Training on {} samples...'.format(len(test_loader.dataset)))
        device = torch.device(cuda_name if torch.cuda.is_available() else "cpu")
        model = modeling().to(device)
        loss_fn = nn.BCELoss()
        optimizer = torch.optim.Adam(model.parameters(), lr=LR)
        best_mse = 1000
        best_ci = 0
        best_epoch = -1
        model_file_name = 'model_' + model_st + '_' + datasetxxxx +  '.model'
        result_file_name = 'result_' + model_st + '_' + datasetxxxx +  '.csv'
        
        model_file_name = 'model_out1.model'
        #model.load_state_dict(torch.load(model_file_name))
        model = torch.load('full_model_out2000.model')
        G,P,N = predicting(model, device, test_loader,TRAIN_BATCH_SIZE)
        auc = aucJ(G, P)
        for i in range(len(G)):
            print (G[i],P[i],N[i])
        print ("AUC:",auc)
        threshold=0.5
        P[P > threshold] = 1
        P[P <= threshold] = 0
        pos_index = P == 1
        pos_index = pos_index.flatten()
        TPR = np.sum(G[pos_index] == 1) * 1.0 / np.sum(pos_index)
        print("TPR: ", TPR)
        precision = np.sum(G[pos_index] == 1) * 1.0 / np.sum(P)
        print("precision: ", precision)
        #auc = aucJ(label, pred)
        accuracy = acc(G, P)
        print ("accuracy: ", accuracy)
        MCC=matthews_corrcoef(G, P)
        print("MCC: ", MCC)

        print ("pos_num: ", pos_num)
        print ("neg_num: ", neg_num)



