import pickle
import numpy as np
# import pysam
from Bio import SeqIO
from multiprocessing import Pool, cpu_count
import time

import random
from typing import List, Union


def save(filename: str, data):
    """
    Saves data to a file using pickle.

    Args:
        filename (str): The name of the file to save to.
        data: The data to be saved.
    """
    with open(filename, "ab") as f:
        pickle.dump(data, f)
        print("Saving complete")


def load(filename):
    """
    Loads data from a file using pickle.

    Args:
        filename (str): The name of the file to load from.

    Returns:
        The loaded data.
    """
    with open(filename, "rb") as f:
        return pickle.load(f)


def sample_multi_sequences(dna_sequence: str, num_samples: int = 10, sample_length: int = 10000) -> Union[
    List[str], str]:
    """
    Samples multiple sequences of a specified length from evenly distributed positions
    within a long DNA sequence and returns them as a list.

    Args:
        dna_sequence (str): The original long DNA sequence.
        num_samples (int): The number of sequences to sample, defaults to 10.
        sample_length (int): The length of each sampled sequence, defaults to 10000.

    Returns:
        Union[List[str], str]: A list containing the sampled DNA sequence strings.
                               If the original sequence is too short, an error message string is returned.
    """
    original_length = len(dna_sequence)

    # Calculate the step size for each sampling interval.
    # Note: this interval is the distance between the start points of the samples.
    interval = original_length // num_samples

    # --- Critical Safety Check ---
    # Check if the original sequence is long enough to ensure the last sampling task
    # does not exceed the end of the sequence.
    # Calculate the starting position of the last sample fragment.
    last_sample_start_pos = (num_samples - 1) * interval
    # If the end position of the last fragment exceeds the length of the original sequence, sampling cannot be completed.
    if last_sample_start_pos + sample_length > original_length:
        required_length = last_sample_start_pos + sample_length
        return (f"Error: The original DNA sequence is too short to complete sampling.\n"
                f"For this configuration, the minimum required sequence length is approximately {required_length} bp, "
                f"but the current sequence length is {original_length} bp.")

    # A list to store all the sampled sequences.
    sampled_sequences_list = []

    # Sample from num_samples evenly distributed positions.
    for i in range(num_samples):
        # Calculate the starting position of the current sample fragment.
        start_pos = i * interval
        # Extract a fragment of the specified length.
        end_pos = start_pos + sample_length
        fragment = dna_sequence[start_pos:end_pos]
        # Add the extracted fragment to the list.
        sampled_sequences_list.append(fragment)

    # Return the list containing all sampled sequences.
    return sampled_sequences_list


def DataCollection(Fasta_Dir, num_samples, sample_length):
    """
    Collects DNA fragments by sampling from a FASTA file.

    Args:
        Fasta_Dir (str): The directory path to the FASTA file.
        num_samples (int): The number of samples to take from each sequence.
        sample_length (int): The length of each sample.

    Returns:
        list: A list of lists, where each inner list contains DNA fragments
              sampled from a sequence in the FASTA file.
    """
    FastaFile = SeqIO.parse(f"{Fasta_Dir}", "fasta")
    DNA_fragments = []
    for Seq in FastaFile:
        seq = str(Seq.seq)
        fragments = sample_multi_sequences(seq, num_samples, sample_length)
        # Ensure that fragments are a list before appending
        if isinstance(fragments, list):
            DNA_fragments.append(fragments)
        else:
            # Print the error message if sampling failed for a sequence
            print(f"Could not sample from sequence {Seq.id}: {fragments}")

    return DNA_fragments



def main(input_fasta):
    data=DataCollection(input_fasta, 10, 10000)
    print(len(data))
    save("SeqFragments.pkl", data)
