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


import time

if __name__=="__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 Step1_Gibbs_and_modification.py <input_file>")
        sys.exit(1)
    file_path = sys.argv[1]
    filename = file_path.split(".")[0]

    with open(file_path, "r") as file:
        dna_sequences = file.read().splitlines()


    output = open(f"output_for_diff_k_{filename}.csv", "w")
    output.write("k,Gibs Sampler Motif, Modified Gibs Sampler Motif, Gibs Sampler Score,Modified Gibs Sampler Score, Gibs Sampler Average Score, Modified Gibs Sampler Average Score, Gibs Sampler Time, Modified Gibs Sampler Time\n")
    gibs_best_score = float('inf')
    gibs_best_motifs = []
    gibs_best_k=0
    modified_gibs_best_score = float('inf')
    modified_gibs_best_motifs = []
    modified_gibs_best_k=0
    t=len(dna_sequences)
    n_iterations = 100

    for k in range(10, 50):  #len(dna_sequences[0]) - 1
        print("k: ", k)
        start_time=time.perf_counter()
        motifs = gibbs_sampler(dna_sequences, k, t, n_iterations)
        stop_time=time.perf_counter()
        gibs_sampler_motif= get_consensus_motif(motifs)
        gibs_sampler_score=compute_score(motifs)
        gibs_sampler_time=stop_time-start_time
        gibs_sampler_avg_score=gibs_sampler_score/k
        print("Gibs Sampler Motifs:")
        for motif in motifs:
            print(motif)
        print("Score: ", gibs_sampler_score)
        print("Average Score: ", gibs_sampler_avg_score)
        if gibs_sampler_avg_score <= gibs_best_score: #if tie, = selects the longest one
            gibs_best_score = gibs_sampler_avg_score
            gibs_best_motifs = motifs
            gibs_best_k=k
        start_time=time.perf_counter()
        motifs = modified_gibbs_sampler(dna_sequences, k, t, n_iterations)
        stop_time=time.perf_counter()
        modified_gibs_sampler_motif= get_consensus_motif(motifs)
        modified_gibs_sampler_score=compute_score(motifs)
        modified_gibs_sampler_time=stop_time-start_time
        modified_gibs_sampler_avg_score=modified_gibs_sampler_score/k
        print("Modified Gibs Sampler Motifs: ")
        for motif in motifs:
            print(motif)
        print("Score: ", modified_gibs_sampler_score)
        print("Average Score: ", modified_gibs_sampler_avg_score)
        if modified_gibs_sampler_avg_score <= modified_gibs_best_score: #if tie, = selects the longest one
            modified_gibs_best_score = modified_gibs_sampler_avg_score
            modified_gibs_best_motifs = motifs
            modified_gibs_best_k=k
        output.write(str(k)+","+gibs_sampler_motif+","+modified_gibs_sampler_motif+","+"{:.4f}".format(gibs_sampler_score)+","+"{:.4f}".format(modified_gibs_sampler_score)+","+"{:.4f}".format(gibs_sampler_avg_score)+","+"{:.4f}".format(modified_gibs_sampler_avg_score)+","+"{:.4f}".format(gibs_sampler_time)+","+"{:.4f}".format(modified_gibs_sampler_time)+"\n")
    print("Gibs Sampler Best Motifs: ")      
    for motif in gibs_best_motifs:
        print(motif)
    print("Score: ")
    print(gibs_best_score)
    print("Modified Gibs Sampler Best Motifs: ")
    for motif in modified_gibs_best_motifs:
        print(motif)
        
    print("Score: ")
    print(modified_gibs_best_score)
    output.close()