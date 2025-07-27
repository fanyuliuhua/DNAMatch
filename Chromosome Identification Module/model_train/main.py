import numpy as np
import torch
import torch.utils.data as Data
from torch.utils.data import DataLoader
import pickle
from Model import LinearNet
import time
# from thop import profile
from parser1 import args
import warnings
import os
import torch.optim as optim
from utils import accuracy
import wandb

# os.environ["WANDB_API_KEY"]=
# os.environ["WANDB_MODE"]="offline"
warnings.filterwarnings('ignore')
# wandb.init(project='LNet')

def load_data(filename):
    with open(filename,"rb") as f:
        return pickle.load(f)

def main():
    score=0
    model=LinearNet()
    model=model.to(args.device)
    optimizer = optim.Adam(model.parameters(),
                           lr=args.lr)
    criterion=torch.nn.CrossEntropyLoss()
    scaler = torch.cuda.amp.GradScaler()
    autocast = torch.cuda.amp.autocast

    data=load_data("data_all.pkl")
    Label=load_data("label_all.pkl")
    for epoch in range(10000):
        for index in range(23):

            # TrainData=load_data(f"dataEmb{index}.pkl")
            # TrainLabel=load_data(f"Label{index}.pkl")
            TrainData=data[index]
            TrainLabel=Label[index]

            # TrainData=np.array([x.cpu().numpy() for x in TrainData])
            # TrainLabel=np.array(TrainLabel)

            print(TrainData.shape)
            print(TrainLabel.shape)

            trainData = torch.from_numpy(TrainData).type(dtype=torch.float32)
            trainData=trainData.to("cuda")
            trainLabel = np.eye(24)[TrainLabel]
            trainLabel = torch.from_numpy(trainLabel).type(dtype=torch.float32)

            torchDataset = Data.TensorDataset(trainData, trainLabel)

            dataset_loader = Data.DataLoader(torchDataset, batch_size=256)
            if index !=22:
                model.train()
                for step ,(x,y) in enumerate(dataset_loader):

                    x=x.to("cuda")
                    y=y.to("cuda")
                    optimizer.zero_grad()
                    with autocast():
                        output=model(x)
                        loss=criterion(output,y)
                    # if epoch%200==0:
                    #
                    #     print(loss,acc_train)
                    # print("loss",loss)
                    scaler.scale(loss).backward()
                    scaler.step(optimizer)
                    scaler.update()
                    acc_train=accuracy(output.float(),y)
                    if step%60==0:
                        wandb.log(
                            {
                                "loss":loss,
                                "acc_train":acc_train
                            }
                        )
            else:
                model.eval()
                with torch.no_grad():
                    acc = []
                    acc2 = []
                    acc3=[]
                    acc4=[]
                    acc5=[]
                    for step, (x, y) in enumerate(dataset_loader):
                        x = x.to("cuda")
                        y = y.to("cuda")
                        output = model(x)
                        acc_val = accuracy(output, y).cpu().detach().numpy()
                        acc2_val = accuracy(output, y,2).cpu().detach().numpy()
                        acc3_val = accuracy(output, y, 3).cpu().detach().numpy()
                        acc4_val = accuracy(output, y, 4).cpu().detach().numpy()
                        acc5_val = accuracy(output, y, 5).cpu().detach().numpy()

                        acc.append(acc_val)
                        acc2.append(acc2_val)
                        acc3.append(acc3_val)
                        acc4.append(acc4_val)
                        acc5.append(acc5_val)
                    wandb.log({
                        "acc_val": sum(acc) / len(acc),
                        "top_2": sum(acc2) / len(acc2),
                        "top_3": sum(acc3) / len(acc3),
                        "top_4": sum(acc4) / len(acc4),
                        "top_5": sum(acc5) / len(acc5),
                    })
                    tmp_score=sum(acc) / len(acc)
                    if score<tmp_score:
                        score=tmp_score
                        torch.save(model.state_dict(),"DNAMatch.pt")
main()



