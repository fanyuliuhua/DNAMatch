from sklearn.metrics import accuracy_score, precision_recall_fscore_support
import pickle
import pysam
import numpy as np


def save(filename:str,data):
    with open(filename,"wb") as f:
        pickle.dump(data,f)
        print("finish")
def load(filename):
    with open(filename,"rb") as f:
        return pickle.load(f)

def getReadMapRate(fileList):
    score=[]
    for file in fileList:
        readfile=pysam.AlignmentFile(f"/Genome/contig/dataBase/result/{file}", "rb", threads=32)
        for read in readfile:
            if not read.query_sequence:
                continue
            if len(read.query_sequence)<1000000:
                continue
            mapped_pairs = read.get_aligned_pairs(matches_only=True)
            seq_len = len(read.query_sequence)
            mapped_len = len(mapped_pairs)
            mapped_rate = mapped_len / seq_len
            score.append(mapped_rate)
    print(len(score))
    print(sum(score)/len(score))
    save("chrAlignmentScore.pkl",score)

def calculate_classification_metrics(y_true: list, y_pred: list, labels=None, average='weighted'):
    """
    Calculates the accuracy, precision, recall, and F1-score for multi-class classification results.

    Args:
    y_true (list): A list of true labels.
    y_pred (list): A list of predicted labels.
    labels (array-like, optional): A list containing all possible class labels. If None,
                                   the labels are inferred from y_true and y_pred.
                                   It is recommended to specify this in multi-class problems
                                   to avoid calculation errors if certain classes are not present in the predictions.
    average (str, optional): The averaging method to use for calculating multi-class metrics. Defaults to 'weighted'.
                             Possible values include:
                             - 'micro': Calculate metrics globally by counting the total true positives,
                                        false negatives, and false positives. Suitable for class imbalance.
                             - 'macro': Calculate metrics for each label, and find their unweighted mean.
                                        Does not take label imbalance into account.
                             - 'weighted': Calculate metrics for each label, and find their average weighted by support
                                           (the number of true instances for each label). Suitable for class imbalance.
                             - 'samples': Calculate metrics for each instance, and find their average
                                          (only for multi-label classification).
                             - None: The scores for each class are returned.

    Returns:
    dict: A dictionary containing the following metrics:
          - 'accuracy': Accuracy score
          - 'precision': Precision score
          - 'recall': Recall score
          - 'f1_score': F1 score
          - 'report': A detailed classification report (if average is not None)
    """
    if not y_true or not y_pred:
        raise ValueError("Input lists y_true and y_pred cannot be empty.")
    if len(y_true) != len(y_pred):
        raise ValueError("The length of y_true and y_pred must be the same.")

    # Ensure inputs are numpy arrays for sklearn functions
    y_true_np = np.array(y_true)
    y_pred_np = np.array(y_pred)

    # 1. Calculate Accuracy
    accuracy = accuracy_score(y_true_np, y_pred_np)

    # 2. Calculate Precision, Recall, and F1-score
    # The precision_recall_fscore_support function returns a tuple containing precision, recall, fscore, support.
    # The returned results will vary depending on the 'average' parameter.
    if labels is None:
        # If labels are not specified, infer all possible classes from the true and predicted values.
        all_labels = np.unique(np.concatenate((y_true_np, y_pred_np)))
        labels_to_use = all_labels
    else:
        labels_to_use = labels

    precision, recall, f1_score, _ = precision_recall_fscore_support(
        y_true_np,
        y_pred_np,
        labels=labels_to_use,
        average=average,
        zero_division=0 # Set the result to 0 when a division by zero occurs.
    )

    # Generate a detailed classification report (optional, but often useful)
    # If 'average' is None, precision, recall, and f1_score are arrays for each class.
    # We need to determine how to present the report based on the 'average' parameter.
    report = None
    if average is None:
        from sklearn.metrics import classification_report
        report = classification_report(y_true_np, y_pred_np, labels=labels_to_use, zero_division=0)

    results = {
        'accuracy': accuracy,
        'precision': precision,
        'recall': recall,
        'f1_score': f1_score,
        'report': report
    }

    return results