import torch
from transformers import AutoTokenizer, AutoModel
import time
import pickle
from torch.utils.data import DataLoader
import warnings
import numpy as np
warnings.filterwarnings('ignore')


def save(filename:str,data):
    with open(filename,"wb") as f:
        pickle.dump(data,f)
        print("存储完成")
def load(filename):
    with open(filename,"rb") as f:
        return pickle.load(f)

def genEmbData(data, tokenizer, model, device):
    model.to(device)
    model.eval()
    seq_embeddings = []
    batch_idx=0
    # # 创建 DataLoader 进行批处理
    # dataloader = DataLoader(data, batch_size=batch_size,num_workers=8)
    start_time = time.time()

    with torch.no_grad():
        for seqs in data:
            seq_embedding=[]
            for seq in seqs:
                batch_idx+=1
                inputs = tokenizer(seq, return_tensors='pt')["input_ids"].cuda()

                with torch.cuda.amp.autocast():
                    hidden_states = model(inputs)[0]
                    # emb = hidden_states.mean(dim=1).cpu().numpy()
                    embedding_mean = torch.mean(hidden_states[0], dim=0)# 对每个序列取平均

                seq_embedding.append(embedding_mean.cpu().numpy())

                if batch_idx % 1000 == 0:
                    elapsed = time.time() - start_time
                    print(f"Processed sequences in {elapsed:.2f} seconds")
            seq_embeddings.append(np.array(seq_embedding))
    # print(len(seq_embeddings))
    # print(len(seq_embeddings[0]))
    # print(len(seq_embeddings[0][0]))
    return seq_embeddings

def main():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    # 预加载 tokenizer 和 model
    tokenizer = AutoTokenizer.from_pretrained(
        "/model/DNABert2",
        trust_remote_code=True
    )
    model = AutoModel.from_pretrained(
        "/model/DNABert2",
        trust_remote_code=True
    )
    data=load("SeqFragments.pkl")
    begintime=time.time()
    embeddings = np.array(genEmbData(data, tokenizer, model, device))

    save(f"SeqFragmentsEmb.pkl",embeddings)
    print("time",time.time()-begintime)


if __name__ == "__main__":
    main()


















