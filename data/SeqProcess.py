import torch
import pickle
import numpy as np
import pysam
import glob


def save(filename:str,data):
    with open(filename,"wb") as f:
        pickle.dump(data,f)
        print("finish")
def load(filename):
    with open(filename,"rb") as f:
        return pickle.load(f)
def create_ans():
    data1=load("outputSV0.pkl")
    data2=load("outputSV1.pkl")
    data3=load("outputSV2.pkl")
    res1=[]
    res2=[]
    res3=[]

    for idx in range(len(data1)):
        res1.extend(data1[idx])
        res2.extend(data2[idx])
        res3.extend(data3[idx])
    ans=[]
    for i in range(len(res1)):
        x1=res1[i]
        x2 = res2[i]
        x3 = res3[i]
        if x1==x2:
            ans.append(x1)
        elif x1==x3:
            ans.append(x1)
        elif x2==x3:
            ans.append(x2)
        else:
            ans.append(-1)
    save("answer.pkl", ans)

def find_most_frequent_number(numbers):
    """
    Returns the number that appears most frequently in a list.

    Args:
        numbers: A list of numbers.

    Returns:
        The number that appears most frequently in the list.
        If the list is empty, returns None.
        If there are multiple numbers with the same highest frequency,
        it returns the first one encountered in the list.
    """
    if not numbers:
        return None

    # Use a dictionary to store the frequency of each number
    frequency = {}
    for number in numbers:
        frequency[number] = frequency.get(number, 0) + 1

    # Find the number with the maximum frequency
    max_frequency = 0
    most_frequent_number = None

    for number, count in frequency.items():
        if count > max_frequency:
            max_frequency = count
            most_frequent_number = number

    return most_frequent_number

class GenFasta():
    def __init__(self, seqs):
        self.seqs = seqs

    def write_fasta(self, output_file, append=False):
        """
        Writes multiple sequences to a file in FASTA format. If append is True, the sequences will be appended to the file.

        Args:
            output_file (str): The name of the output file.
            append (bool): Whether to open the file in append mode. Defaults to False (overwrite mode).
        """
        mode = 'a' if append else 'w'
        with open(output_file, mode) as fasta_file:
            for identifier, sequence in self.seqs:
                print(identifier, len(sequence))
                # Write the sequence identifier
                fasta_file.write(f">{identifier}\n")
                # Write the sequence in lines with a fixed width (e.g., 60 characters)
                for i in range(0, len(sequence), 60):
                    fasta_file.write(sequence[i:i + 60] + '\n')


def PreProcessData():
    ans = load("answer.pkl")
    dataBase = {
        0: [],
        1: [],
        2: [],
        3: [],
        4: [],
        5: [],
        6: [],
        7: [],
        8: [],
        9: [],
        10: [],
        11: [],
        12: [],
        13: [],
        14: [],
        15: [],
        16: [],
        17: [],
        18: [],
        19: [],
        20: [],
        21: [],
        22: [],
        23: [],
        -1: [],
    }
    from Bio import SeqIO
    bamfile = SeqIO.parse("Seqs.fasta", "fasta")
    idx_loc = 0
    count = 0
    batchIdx=0
    for read in bamfile:
        if idx_loc%10000==0:
            print(idx_loc)
        count+=1
        print(idx_loc)
        cls=ans[idx_loc]
        dataBase[cls].append((idx_loc,str(read.seq.upper())))
        idx_loc+=1
        count+=1
        if count%40000==0:
            for k, v in dataBase.items():
                genData = GenFasta(dataBase[k])
                genData.write_fasta(f"/io/fanyuliuhua/Genome/contig/dataBase/data/@chr{k+1}@_{batchIdx}.fasta")
                dataBase[k]=[]
        batchIdx+=1
    for k,v in dataBase.items():
        print(idx_loc)
        genData = GenFasta(dataBase[k])
        genData.write_fasta(f"/io/fanyuliuhua/Genome/contig/dataBase/data/@chr{k+1}@end.fasta")
        dataBase[k] = []



def concat_fasta():
    #42937 42947
    fileName=["WGS_SV_Million_h1.sam","WGS_SV_Million_h2.sam"]
    idx_loc = 0
    data=[]
    label=[]
    for file in fileName:
        bamfile = pysam.AlignmentFile(f"/home/fanyuliuhua/github/SV/{file}", "rb", threads=24)
        for read in bamfile:
            data.append((idx_loc,read.query_sequence))
            label.append(read.reference_name)
            idx_loc+=1
            if len(data)==10000:
                print(idx_loc)
                gendata=GenFasta(data)
                gendata.write_fasta("/io/fanyuliuhua/SV/bamfile.fasta",append=True)
                data=[]
        if len(data)>0:
            gendata = GenFasta(data)
            gendata.write_fasta("/io/fanyuliuhua/SV/bamfile.fasta", append=True)
            data = []
    save("label.pkl", label)
    print("finish",idx_loc)



# PreProcessData()
def eval():
    # from collections import Counter
    ans = load("answer.pkl")
    label=load("label.pkl")
    # count=Counter(label)
    # for k,v in count.items():
    #     print(k,v)
    check=[]
    for i in range(len(ans)):
        if ans[i]==label[i]:
            check.append(1)
        elif ans[i]==-1:
            check.append(1)
        else:
            check.append(0)
    print(sum(check)/len(check))
def concat():
    chrName=["chr0","chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14"
             ,"chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chr23","chr24"]
    for chr in chrName:
        files=glob.glob(f"/Genome/contig/dataBase/data/@{chr}@*.fasta")
        output_file = f"{chr}.fasta"
        # print(files)
        with open("/Genome/contig/dataBase/SeqData/"+output_file, "w") as outfile:
            for fname in files:
                with open(fname, "r") as infile:
                    outfile.write(infile.read())


# def main():
#     eval()
#     create_ans()
#     PreProcessData()
#     concat()
#     concat_fasta()
# main()




