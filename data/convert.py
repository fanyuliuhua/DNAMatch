import torch
from transformers import AutoTokenizer, AutoModel
import time
import pickle
import numpy as np
from torch.utils.data import DataLoader
import warnings
import collections


warnings.filterwarnings('ignore')

def save(filename:str,data):
    with open(filename,"wb") as f:
        pickle.dump(data,f)
        print("存储完成")
def load(filename):
    with open(filename,"rb") as f:
        return pickle.load(f)

file_dir="D:\\PyProjects\\data\\"
strToNum={
    "chr1":0,
    "chr2":1,
    "chr3":2,
    "chr4":3,
    "chr5": 4,
    "chr6": 5,
    "chr7": 6,
    "chr8": 7,
    "chr9": 8,
    "chr10": 9,
    "chr11": 10,
    "chr12": 11,
    "chr13": 12,
    "chr14": 13,
    "chr15": 14,
    "chr16": 15,
    "chr17": 16,
    "chr18": 17,
    "chr19": 18,
    "chr20": 19,
    "chr21": 20,
    "chr22": 21,
    "chrX": 22,
    "chrY": 23,
}


#GUOYING
strToNum={
    "NC_004353.4":0,
    "NC_004354.4":1,
    "NC_024512.1":2,
    "NT_033777.3":3,
    "NT_033778.4":4,
    "NT_033779.5":5,
    "NT_037436.4":6,
}



strToNum={
    "NC_003070.9_Arabidopsis_thaliana_chromosome_1_sequence":0,
    "NC_003076.8_Arabidopsis_thaliana_chromosome_5,_partial_sequence":1,
    "NC_003074.8_Arabidopsis_thaliana_chromosome_3,_partial_sequence":2,
    "NC_003071.7_Arabidopsis_thaliana_chromosome_2,_partial_sequence":3,
    "NC_003075.7_Arabidopsis_thaliana_chromosome_4,_partial_sequence":4,
}


def mapDict(x):
    return strToNum[x]

# x=load("label.pkl")
# print(x[:10])
# document=f"./Genome/ninanjie/label.pkl"
# data=load(document)
# result=list(map(mapDict,data))
# result=np.array(result)
# save(f"./Genome/ninanjie/labelN.pkl",result)
# for i in range(10):
#     x=load(f"./dataBase/output{i}.pkl")
#     print(x[:10])
# y=load(f"./dataBase/label.pkl")
# print(y[:10])