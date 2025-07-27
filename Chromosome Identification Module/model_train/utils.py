import numpy as np
import torch

def accuracy(output, labels,top_k=1):
    pred=torch.topk(output,top_k,dim=1).indices
    val=torch.argmax(labels,1)
    correct=pred.eq(val.view(-1,1))
    correct = correct.sum(dim=1)
    return correct.float().mean()
