import random
import math
import numpy as np
import sys

def make_Profile(motifs , k , t , stringList):
    count = [[],[],[],[]]
    for i in range(k):
        count[0].append(1)
        count[1].append(1)
        count[2].append(1)
        count[3].append(1)
    
    for column in range(k):
        for row in range(t):
            
            temp = motifs[row] + column
            
            if stringList[row][temp] == 'A':
                count[0][column] +=1
            elif stringList[row][temp] == 'C':
                count[1][column] +=1
            elif stringList[row][temp] == 'G':
                count[2][column] +=1
            else:
                count[3][column] +=1
    
    return count


def make_motifs(profile , k , t ,stringList, stringLength):
    motifs = []
    for row in range(t):
        max = -1
        num_index = -1
        for column in range(stringLength-k+1):
            multiplication=1
            
            for i in range(k):
                temp = stringList[row][column + i]
                if temp == 'A' :
                    multiplication = multiplication * profile[0][i]
                elif temp == 'C' :
                    multiplication = multiplication * profile[1][i]
                elif temp == 'G' :
                    multiplication = multiplication * profile[2][i]
                else:
                    multiplication = multiplication * profile[3][i]

            if(max < multiplication):
                max = multiplication
                num_index = column
        motifs.append(num_index)
    
    return motifs
            
def make_motifs_modification(profile , k , t ,stringList, stringLength):
    motifs = []
    for row in range(t):
        
        num_index = -1
        vector=[]
        for column in range(stringLength-k+1):
            multiplication=1
            list1 = []
            for i in range(k):
                temp = stringList[row][column + i]
                if temp == 'A' :
                    multiplication = multiplication * profile[0][i]
                elif temp == 'C' :
                    multiplication = multiplication * profile[1][i]
                elif temp == 'G' :
                    multiplication = multiplication * profile[2][i]
                else:
                    multiplication = multiplication * profile[3][i]

            list1.append(multiplication)
            list1.append(column)
            vector.append(list1)

            
        
        vector.sort(reverse=True)
        random_number = random.randint(0,3)
        num_index = vector[random_number][1]

        motifs.append(num_index)
    
    return motifs

def score_function_by_entropy(motifs,k,t,stringList):
     profile = make_Profile(motifs,k,t,stringList)
     entropy=0.0
     for column in range(k):
         sum=0
         for row in range(4):
             sum +=profile[row][column]
         for row in range(4):
             temp = profile[row][column] / sum
             temp = -temp*np.log2(temp)
             entropy +=temp
     return entropy
          

def get_consensus_motif(motifs, k,t , stringList):
    profile = make_Profile(motifs,k,t,stringList)
    consensus = ""
    for column in range(k):
        max = -1
        
        for row in range(4):
            if profile[row][column] > max:
                max = profile[row][column]
        
        for row in range(4):
            if profile[row][column] == max:
                if row == 0:
                    consensus +="A"
                elif row == 1:
                    consensus +="C"
                elif row == 2:
                    consensus +="G"
                else:
                    consensus +="T"
             
                break

    return consensus


def score_function_by_hamming_distance(motifs, k,t , stringList):
    profile = make_Profile(motifs,k,t,stringList)
    sum=0
    for column in range(k):
        max = -1
        
        for row in range(4):
            if profile[row][column] > max:
                max = profile[row][column]
        
        for row in range(4):
            if profile[row][column] != max:
                sum += profile[row][column]
                sum -=1
    return sum

def local_optimization_for_entropy(motifs, k, t, stringList, stringLength):
    best_motifs = motifs[:]
    best_score = score_function_by_entropy(best_motifs, k, t, stringList)
    improved = False
    
    while improved == False:
        
        
        for i in range(t):
            original_motif = motifs[i]
            count = 0
            for j in range(stringLength - k + 1):
                motifs[i] = j
                score = score_function_by_entropy(motifs, k, t, stringList)
                if score < best_score:
                    best_motifs = motifs[:]
                    best_score = score
                    improved = True
                    count +=1
                    if count == 2:
                        break 
            motifs[i] = original_motif
            if improved == True:
                break
               
            
    
    return best_motifs

def local_optimization_for_hamming(motifs, k, t, stringList, stringLength):
    best_motifs = motifs[:]
    best_score = score_function_by_hamming_distance(best_motifs, k, t, stringList)
    improved = False
    
    while improved == False:
        
        
        for i in range(t):
            original_motif = motifs[i]
            count = 0
            for j in range(stringLength - k + 1):
                motifs[i] = j
                score = score_function_by_hamming_distance(motifs, k, t, stringList)
                if score < best_score:
                    best_motifs = motifs[:]
                    best_score = score
                    improved = True
                    count +=1
                    if count == 2:
                        break 
            motifs[i] = original_motif
            if improved == True:
                break
               
            
    
    return best_motifs



def Randomized_Motif_Search(stringList, k , t , stringLength):
    motifs = []
    for i in range(t):
        random_number = random.randint(0,stringLength-k)
        motifs.append(random_number)
    
    bestMotif = motifs
    iter = 0

    while True:
        iter +=1
        profile = make_Profile(motifs,k,t , stringList)

        motifs = make_motifs(profile,k,t,stringList,stringLength)
        

        if(score_function_by_entropy(motifs,k,t,stringList) < score_function_by_entropy(bestMotif,k,t,stringList)):
            bestMotif = motifs
        else:
            # print("Iteration for this call : ",iter)
            return bestMotif

def Randomized_Motif_Search_With_hamming(stringList, k , t , stringLength):
    motifs = []
    for i in range(t):
        random_number = random.randint(0,stringLength-k)
        motifs.append(random_number)
    
    bestMotif = motifs
    iter = 0

    while True:
        iter +=1
        profile = make_Profile(motifs,k,t , stringList)

        motifs = make_motifs(profile,k,t,stringList,stringLength)
        

        if(score_function_by_hamming_distance(motifs,k,t,stringList) < score_function_by_hamming_distance(bestMotif,k,t,stringList)):
            bestMotif = motifs
        else:
            # print("Iteration for this call : ",iter)
            return bestMotif


def Modified_Randomized_Motif_Search(stringList, k , t , stringLength):
    motifs = []
    for i in range(t):
        random_number = random.randint(0,stringLength-k)
        motifs.append(random_number)
    
    bestMotif = motifs
    iter = 0

    while True:
        iter +=1
        profile = make_Profile(motifs,k,t , stringList)

        motifs = make_motifs_modification(profile,k,t,stringList,stringLength)
        motifs = local_optimization_for_entropy(motifs, k, t, stringList, stringLength)

        if(score_function_by_entropy(motifs,k,t,stringList) < score_function_by_entropy(bestMotif,k,t,stringList)):
            bestMotif = motifs
        else:
            # print("Iteration for this call : ",iter)
            return bestMotif

def Modified_Randomized_Motif_Search_With_Hamming(stringList, k , t , stringLength):
    motifs = []
    for i in range(t):
        random_number = random.randint(0,stringLength-k)
        motifs.append(random_number)
    
    bestMotif = motifs
    iter = 0

    while True:
        iter +=1
        profile = make_Profile(motifs,k,t , stringList)

        motifs = make_motifs_modification(profile,k,t,stringList,stringLength)
        motifs = local_optimization_for_hamming(motifs, k, t, stringList, stringLength)

        if(score_function_by_hamming_distance(motifs,k,t,stringList) < score_function_by_hamming_distance(bestMotif,k,t,stringList)):
            bestMotif = motifs
        else:
            # print("Iteration for this call : ",iter)
            return bestMotif

def run_many_times(function_number,iteration,stringList,k,totalString,stringLength):
    best_score = -1
    best_motif=[]
    
    for i in range(iteration):
        resultMotif = []
        if function_number == 1:
            resultMotif =Randomized_Motif_Search(stringList,k,totalString,stringLength)
        elif function_number == 2:
            resultMotif = Randomized_Motif_Search_With_hamming(stringList, k , t , stringLength)
        elif function_number == 3:
            resultMotif = Modified_Randomized_Motif_Search(stringList, k , t , stringLength)
        else:
            resultMotif =  Modified_Randomized_Motif_Search_With_Hamming(stringList, k , t , stringLength)

        if i == 0:
            best_score=score_function_by_entropy(resultMotif,k,totalString,stringList)
            best_motif = resultMotif
            
            
        else:
            tempScore = score_function_by_entropy(resultMotif,k,totalString,stringList)
            if best_score > tempScore:
                best_score = tempScore
                best_motif = resultMotif

    return best_motif            
 

def comparison(stringList,k,t,stringLength,output):
    motifs = run_many_times(1,3,stringList,k,t,stringLength)
    Randomized_motif_with_entropy = get_consensus_motif(motifs,k,t,stringList)
    Randomized_score_with_entropy = score_function_by_entropy(motifs,k,t,stringList)
    motifs=run_many_times(2,3,stringList,k,t,stringLength)
    Randomized_motif_with_hamming = get_consensus_motif(motifs,k,t,stringList)
    Randomized_score_with_hamming = score_function_by_entropy(motifs,k,t,stringList)
    motifs = run_many_times(3,3,stringList,k,t,stringLength)
    Modified_Randomized_motif_with_entropy = get_consensus_motif(motifs,k,t,stringList)
    Modified_Randomized_score_with_entropy = score_function_by_entropy(motifs,k,t,stringList)
    motifs=run_many_times(4,3,stringList,k,t,stringLength)
    Modified_Randomized_motif_with_hamming = get_consensus_motif(motifs,k,t,stringList)
    Modified_Randomized_score_with_hamming = score_function_by_entropy(motifs,k,t,stringList)
    kmer = f'{k}'
    output.write(kmer+","+Randomized_motif_with_entropy+","+Randomized_motif_with_hamming+","+Modified_Randomized_motif_with_entropy+","+Modified_Randomized_motif_with_hamming+","+"{:.4f}".format(Randomized_score_with_entropy)+","+"{:.4f}".format(Modified_Randomized_score_with_entropy)+","+"{:.4f}".format(Modified_Randomized_score_with_hamming)+","+"{:.4f}".format(Modified_Randomized_score_with_hamming)+"\n")


if __name__=="__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 Step2_Randomized_and_modification.py <input_file> 1st_k_value 2nd_k_value")
        sys.exit(1)
    file_path = sys.argv[1]
    filename = file_path.split(".")[0]
    f = open(file_path, "r")

    totalString = 0
    stringLength = 0
    stringList = []

    while True:
        str = f.readline()
        
        if len(str) == 0:
            break
        #ignore the newline character
        str = str[:-1]
        stringList.append(str)
        totalString += 1
        stringLength = len(str)

    f.close()
    t=totalString
    
    output=open("output_for_randomized_comparison_"+file_path[:-3]+"csv", "w")
    output.write("k,Randomized Motif With Entropy, Randomized Motif With Hamming, Modified Randomized Motif With Entropy, Modified Randomized Motif With Hamming, Randomized Score With Entropy,Modified Randomized Score With Entropy,  Randomized Score With Hamming,Modified Randomized Score With Hamming\n")
    
    comparison(stringList,int(sys.argv[2]),t,stringLength,output)
    comparison(stringList,int(sys.argv[3]),t,stringLength,output)
    output.close

