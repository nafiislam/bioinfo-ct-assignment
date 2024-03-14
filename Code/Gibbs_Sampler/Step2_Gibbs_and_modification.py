import numpy as np
import random
import sys

def compute_profile_matrix(motifs):
    profile_matrix = []
    for i in range(len(motifs[0])):
        column = [motif[i] for motif in motifs]
        profile_matrix.append({
            'A': (column.count('A')+1) / (len(column)+4),
            'C': (column.count('C')+1) / (len(column)+4),
            'G': (column.count('G')+1) / (len(column)+4),
            'T': (column.count('T')+1) / (len(column)+4)
        })
    return profile_matrix

# compute the score from Entropy
def compute_score(motifs):
    profile_matrix = compute_profile_matrix(motifs)
    entropy = []
    
    # print(len(profile_matrix))
    
    for i in range(len(profile_matrix)):
        if(profile_matrix[i]['A'] != 0):
            entropy.append(-profile_matrix[i]['A'] * np.log2(profile_matrix[i]['A']))
        if(profile_matrix[i]['C'] != 0):
            entropy.append(-profile_matrix[i]['C'] * np.log2(profile_matrix[i]['C']))
        if(profile_matrix[i]['G'] != 0):
            entropy.append(-profile_matrix[i]['G'] * np.log2(profile_matrix[i]['G']))
        if(profile_matrix[i]['T'] != 0):
            entropy.append(-profile_matrix[i]['T'] * np.log2(profile_matrix[i]['T']))            
        
    return sum(entropy)

def compute_score1(motifs):
    score=0
    consensus=get_consensus_motif(motifs)
    for motif in motifs:
        score+=get_hamming_distance(consensus, motif)
    return score

def calculate_prob_kmer(dna_sequence, k, profile_matrix):
    probabilities = []
    for i in range(len(dna_sequence) - k + 1):
        kmer = dna_sequence[i:i + k]
        probability = 1
        for j in range(k):
            probability *= profile_matrix[j][kmer[j]]
        probabilities.append(probability)
    return probabilities


def select_substring_from_probablities(dna_sequence, probabilities, k):
    min_value = float('inf')  # Initialize min_value with positive infinity
    for i in range(len(probabilities)):
        if(probabilities[i] != 0 and probabilities[i] < min_value):
            min_value = probabilities[i]
    
    for i in range(len(probabilities)):
        probabilities[i] = probabilities[i] / min_value
        
    # normalize the probablities which are not zero
    sum_probablities = sum(probabilities)
    if(sum_probablities != 0):
        for i in range(len(probabilities)):
            if(probabilities[i] != 0):
                probabilities[i] = probabilities[i] / sum_probablities
                
        
        sorted_probabilities = []
        sorted_probabilities_index = []
        temp_probabilities = probabilities
        
        for i in range(len(temp_probabilities)):
            min_val = float('inf')
            min_index = -1
            for j in range(len(temp_probabilities)):
                if(temp_probabilities[j] < min_val):
                    min_val = temp_probabilities[j]
                    min_index = j
            sorted_probabilities.append(min_val)
            sorted_probabilities_index.append(min_index)
            temp_probabilities[min_index] = float('inf')
            
            
        # print(sorted_probabilities_index)
        
        # print(min(sorted_probabilities_index))
        # print(max(sorted_probabilities_index))
        
        # print(sum(sorted_probabilities))
                
        random_number = random.uniform(min(sorted_probabilities), max(sorted_probabilities))
        for i in range(len(sorted_probabilities)):
            if(random_number < sorted_probabilities[i]):
                return dna_sequence[sorted_probabilities_index[i]:sorted_probabilities_index[i] + k]
            
    else:
        random_index = random.randint(0, len(dna_sequence) - k)
        return dna_sequence[random_index:random_index + k]
        
        

def get_consensus_motif(motifs):
    consensus = ""
    for i in range(len(motifs[0])):
        column = [motif[i] for motif in motifs]
        consensus += max(set(column), key=column.count)
    return consensus

def get_hamming_distance(motif1, motif2):
    distance=0
    for i in range(len(motif1)):
        if motif1[i] != motif2[i]:
            distance+=1
    return distance

def get_worst_k_mer(motifs):
    hamming_distances=[]
    worst_k_mers=[]
    worst_distances=0
    for i in range(len(motifs)):
        temp_motifs=motifs[:i]+motifs[i+1:]
        hamming_distances.append(get_hamming_distance(get_consensus_motif(temp_motifs), motifs[i]))
        if hamming_distances[-1]>worst_distances:
            worst_distances=hamming_distances[-1]
    for i in range(len(hamming_distances)):
        if hamming_distances[i]==worst_distances:
            worst_k_mers.append(i)
    return worst_k_mers[random.randint(0,len(worst_k_mers)-1)]

def gibbs_sampler(dna_sequences, k, t, n_iterations):
    motifs = []
    best_motifs = []
    # randomly select starting positions from every sequence
    for i in range(len(dna_sequences)):
        start = random.randint(0, len(dna_sequences[i]) - k)
        motifs.append(dna_sequences[i][start:start + k])
    best_motifs = motifs
    
    for j in range(n_iterations):
        i = random.randint(0, t - 1)
        motifs.pop(i)
        profile_matrix = compute_profile_matrix(motifs)
        probabilities = calculate_prob_kmer(dna_sequences[i], k, profile_matrix)
        # print(len(probabilities))
        selected_substring = select_substring_from_probablities(dna_sequences[i], probabilities, k)
        motifs = motifs[:i] + [selected_substring] + motifs[i:]
        
        # print(selected_substring)
                
        # for motif in motifs:
        #     print(motif)
        
        if(compute_score(motifs) < compute_score(best_motifs)):
            best_motifs = motifs
            
    return best_motifs   

def gibbs_sampler_with_hamming(dna_sequences, k, t, n_iterations):
    motifs = []
    best_motifs = []
    # randomly select starting positions from every sequence
    for i in range(len(dna_sequences)):
        start = random.randint(0, len(dna_sequences[i]) - k)
        motifs.append(dna_sequences[i][start:start + k])
    best_motifs = motifs
    
    for j in range(n_iterations):
        i = random.randint(0, t - 1)
        motifs.pop(i)
        profile_matrix = compute_profile_matrix(motifs)
        probabilities = calculate_prob_kmer(dna_sequences[i], k, profile_matrix)
        # print(len(probabilities))
        selected_substring = select_substring_from_probablities(dna_sequences[i], probabilities, k)
        motifs = motifs[:i] + [selected_substring] + motifs[i:]
        
        # print(selected_substring)
                
        # for motif in motifs:
        #     print(motif)
        
        if(compute_score1(motifs) < compute_score1(best_motifs)):
            best_motifs = motifs
            
    return best_motifs   

def modified_gibbs_sampler(dna_sequences, k, t, n_iterations):
    motifs = []
    best_motifs = []
    # randomly select starting positions from every sequence
    for i in range(len(dna_sequences)):
        start = random.randint(0, len(dna_sequences[i]) - k)
        motifs.append(dna_sequences[i][start:start + k])
    best_motifs = motifs
    
    for j in range(n_iterations):
        i = get_worst_k_mer(motifs)
        motifs.pop(i)
        profile_matrix = compute_profile_matrix(motifs)
        probabilities = calculate_prob_kmer(dna_sequences[i], k, profile_matrix)
        # print(len(probabilities))
        selected_substring = select_substring_from_probablities(dna_sequences[i], probabilities, k)
        motifs = motifs[:i] + [selected_substring] + motifs[i:]
        
        # print(selected_substring)
                
        # for motif in motifs:
        #     print(motif)
        
        if(compute_score(motifs) < compute_score(best_motifs)):
            best_motifs = motifs
            
    return best_motifs   


def modified_gibbs_sampler_with_hamming(dna_sequences, k, t, n_iterations):
    motifs = []
    best_motifs = []
    # randomly select starting positions from every sequence
    for i in range(len(dna_sequences)):
        start = random.randint(0, len(dna_sequences[i]) - k)
        motifs.append(dna_sequences[i][start:start + k])
    best_motifs = motifs
    
    for j in range(n_iterations):
        i = get_worst_k_mer(motifs)
        motifs.pop(i)
        profile_matrix = compute_profile_matrix(motifs)
        probabilities = calculate_prob_kmer(dna_sequences[i], k, profile_matrix)
        # print(len(probabilities))
        selected_substring = select_substring_from_probablities(dna_sequences[i], probabilities, k)
        motifs = motifs[:i] + [selected_substring] + motifs[i:]
        
        # print(selected_substring)
                
        # for motif in motifs:
        #     print(motif)
        
        if(compute_score1(motifs) < compute_score1(best_motifs)):
            best_motifs = motifs
            
    return best_motifs   

def comparison(dna_sequences, k, t, n_iterations, output):
    motifs=gibbs_sampler(dna_sequences, k, t, n_iterations)
    gibs_sampler_motif_with_entropy=get_consensus_motif(motifs)
    gibs_sampler_score_with_entropy=compute_score(motifs)
    motifs=gibbs_sampler_with_hamming(dna_sequences, k, t, n_iterations)
    gibs_sampler_motif_with_hamming=get_consensus_motif(motifs)
    gibs_sampler_score_with_hamming=compute_score(motifs)
    motifs=modified_gibbs_sampler(dna_sequences, k, t, n_iterations)
    modified_gibs_sampler_motif_with_entropy=get_consensus_motif(motifs)
    modified_gibs_sampler_score_with_entropy=compute_score(motifs)
    motifs=modified_gibbs_sampler_with_hamming(dna_sequences, k, t, n_iterations)
    modified_gibs_sampler_motif_with_hamming=get_consensus_motif(motifs)
    modified_gibs_sampler_score_with_hamming=compute_score(motifs)
    output.write(str(k)+","+gibs_sampler_motif_with_entropy+","+gibs_sampler_motif_with_hamming+","+modified_gibs_sampler_motif_with_entropy+","+modified_gibs_sampler_motif_with_hamming+","+"{:.4f}".format(gibs_sampler_score_with_entropy)+","+"{:.4f}".format(modified_gibs_sampler_score_with_entropy)+","+"{:.4f}".format(gibs_sampler_score_with_hamming)+","+"{:.4f}".format(modified_gibs_sampler_score_with_hamming)+"\n")


import time
import sys

if __name__=="__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 Step2_Gibbs_and_modification.py <input_file> <k_value1> <k_value2>")
        sys.exit(1)
    file_path = sys.argv[1]

    with open(file_path, "r") as file:
        dna_sequences = file.read().splitlines()

    t=len(dna_sequences)
    n_iterations = 100

    output=open("output_for_comparison_"+file_path[:-3]+"csv", "w")
    output.write("k,Gibs Sampler Motif With Entropy, Gibs Sampler Motif With Hamming, Modified Gibs Sampler Motif With Entropy, Modified Gibs Sampler Motif With Hamming, Gibs Sampler Score With Entropy,Modified Gibs Sampler Score With Entropy,  Gibs Sampler Score With Hamming,Modified Gibs Sampler Score With Hamming\n")
    comparison(dna_sequences, int(sys.argv[2]), t, n_iterations, output)
    comparison(dna_sequences, int(sys.argv[3]), t, n_iterations, output)
    output.close()
    