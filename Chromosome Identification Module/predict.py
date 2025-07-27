import numpy as np
import torch
import torch.utils.data as Data
import sys

sys.path.append("/home/fanyuliuhua/DnaMatch")
from torch.utils.data import DataLoader
import pickle
from model_train.Model import LinearNet
import time
# from thop import profile
from model_train.parser1 import args
import warnings
import os
import torch.optim as optim
from model_train.utils import accuracy
import wandb
def save(filename:str,data):
    with open(filename,"wb") as f:
        pickle.dump(data,f)
        print("存储完成")
def load_data(filename):
    with open(filename,"rb") as f:
        return pickle.load(f)
model=LinearNet()
model.load_state_dict(torch.load("DNAMatch.pt"))
model.to("cuda")
dataFiles=["SeqReadsEmb.pkl"]
for idx in range(len(dataFiles)):
    data = load_data(f"{dataFiles[idx]}")
    print(len(data))
    trainData = torch.from_numpy(data).type(dtype=torch.float32)
    trainData = trainData.to("cuda")
    acc=[]
    model.eval()
    with torch.no_grad():
        out=[]
        for x in trainData:
            output = model(x)
            out.append(torch.argmax(output,dim=1).cpu().detach().numpy().tolist())
    save("predict.pkl",out)

