# import vcf
import pickle
import sys
import time


def save(filename:str,data):
    with open(filename,"wb") as f:
        pickle.dump(data,f)
        print("存储完成")
def load(filename):
    with open(filename,"rb") as f:
        return pickle.load(f)
def read_fa_file(file_path):
    with open(file_path, "r") as file:
        next(file)
        sequence = file.read().replace("\n", "")  # 读取整个文件并移除换行符
    return sequence

def get_N_pos():
    idx=[1,2,3,4,5,6,7,8,9,10,
         11,12,13,14,15,16,17,18,19,20,21,22,"X","Y"]
    # idx=[14]
    Npos=dict()
    for id in idx:
        file = f"D:\\data\\refg\\Homo_sapiens.GRCh37.dna.chromosome.{id}.fa"
        seq=read_fa_file(file)
        print(len(seq))
        start=-1
        end=-1
        for i in range(len(seq)):
            if seq[i]=="N":
                if start==-1:
                    start=i
                    end=i
                elif i==len(seq)-1:
                    end+=1
                    Npos.setdefault(id, []).append((start, end))
                else:
                    end+=1
            else:
                if start!=-1:
                    Npos.setdefault(id,[]).append((start,end))
                    start=-1
                    end=-1
                else:
                    continue
        print(Npos[id])#有碱基N的区域
        print(len(Npos[id]))
        # print(seq[148361358:148511357+1])
    save("unknownBase.pkl",Npos)
def filter_UnknownBase():
    idx=[1,2,3,4,5,6,7,8,9,10,
         11,12,13,14,15,16,17,18,19,20,21,22,"X","Y"]
    signal_UnknownBase=load("unknownBase.pkl")
    KnownBase=dict()
    for id in idx:
        file = f"D:\\data\\refg\\Homo_sapiens.GRCh37.dna.chromosome.{id}.fa"
        seq=read_fa_file(file)
        chr_UnknownBase=signal_UnknownBase[id]
        chr_UnknownBase=sorted(chr_UnknownBase,key=lambda x:x[0],reverse=False)
        seqlen=len(seq)
        print(seqlen)
        pre=None
        for num in range(len(chr_UnknownBase)):
            domain=chr_UnknownBase[num]
            if num==0:
                if num==len(chr_UnknownBase)-1:
                    KnownBase.setdefault(id, []).append((0, domain[0]))
                    KnownBase.setdefault(id, []).append((domain[0], seqlen))
                else:
                    KnownBase.setdefault(id, []).append((0, domain[0]))
                    pre = domain[1]
            elif num==len(chr_UnknownBase)-1:
                KnownBase[id].append((pre,domain[0]))
                KnownBase[id].append((domain[1],seqlen))
                print(f"rainflow{id}")
            else:
                KnownBase[id].append((pre,domain[0]))
                pre=domain[1]
    print(KnownBase)
    save("KnownBase.pkl",KnownBase)

def concat_chr():
    count=0
    data=load("KnownBaseFilter.pkl")
    for k,v in data.items():
        init_seq=""
        file = f"D:\\data\\refg\\Homo_sapiens.GRCh37.dna.chromosome.{k}.fa"
        seq = read_fa_file(file)
        for tmp in v:
            m,n=tmp
            if n-m<=1:
                continue
            init_seq+=seq[m+1:n]
        print(len(init_seq))
        count+=len(init_seq)
        if "N" in init_seq:
            print("错误")
        save(f"./seqs/chr{k}.pkl",init_seq)