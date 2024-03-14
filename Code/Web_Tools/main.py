import numpy as np
import random
import sys


def compute_profile_matrix(motifs):
    profile_matrix = []
    for i in range(len(motifs[0])):
        column = [motif[i] for motif in motifs]
        profile_matrix.append({
            'A': (column.count('A') + 1) / (len(column) + 4),
            'C': (column.count('C') + 1) / (len(column) + 4),
            'G': (column.count('G') + 1) / (len(column) + 4),
            'T': (column.count('T') + 1) / (len(column) + 4)
        })
    return profile_matrix


# compute the score from Entropy
def compute_score(motifs):
    profile_matrix = compute_profile_matrix(motifs)
    entropy = []

    # print(len(profile_matrix))

    for i in range(len(profile_matrix)):
        if profile_matrix[i]['A'] != 0:
            entropy.append(-profile_matrix[i]['A'] * np.log2(profile_matrix[i]['A']))
        if profile_matrix[i]['C'] != 0:
            entropy.append(-profile_matrix[i]['C'] * np.log2(profile_matrix[i]['C']))
        if profile_matrix[i]['G'] != 0:
            entropy.append(-profile_matrix[i]['G'] * np.log2(profile_matrix[i]['G']))
        if profile_matrix[i]['T'] != 0:
            entropy.append(-profile_matrix[i]['T'] * np.log2(profile_matrix[i]['T']))

    return sum(entropy)


def get_consensus_motif(motifs):
    consensus = ""
    for i in range(len(motifs[0])):
        column = [motif[i] for motif in motifs]
        consensus += max(set(column), key=column.count)
    return consensus


def get_hamming_distance(motif1, motif2):
    distance = 0
    for i in range(len(motif1)):
        if motif1[i] != motif2[i]:
            distance += 1
    return distance


def compute_score1(motifs):
    score = 0
    consensus = get_consensus_motif(motifs)
    for motif in motifs:
        score += get_hamming_distance(consensus, motif)
    return score




allFiles = ['yst04r_MEME.txt', 'yst04r_MEMEChIP.txt', 'yst08r_MEME.txt', 'yst08r_MEMEChIP.txt', 'hm03_MEME.txt', 'hm03_MEMEChIP.txt',]

csv = open('output.csv', 'w')
csv.write("k,Data, Tool,Motif,Scoring function, Score, Avg Score\n")
for file in allFiles:
    with open(file, 'r') as f:
        data = file.split('_')[0]
        print(data)
        tool = file.split('_')[1].split('.')[0]
        print(tool)

        lines = f.readlines()
        motifSet = []
        for line in lines:
            motifSet.append(line.strip())
        entropy_score = compute_score(motifSet)
        avg_entropy_score = entropy_score / len(motifSet[0])

        csv.write(str(len(motifSet[0])) + "," + data + "," + tool + "," + get_consensus_motif(motifSet) +",Entropy," + "{:.4f}".format(entropy_score) + ",{:.4f}".format(avg_entropy_score) + "\n")


