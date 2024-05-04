import os
import numpy as np

#cwd = os.path.dirname(os.path.abspath("C:\Users\Vonag\Downloads\Project\Project\datapoints\in1.txt"))
cwd = r"C:\Users\Vonag\Downloads\Project\Project"
print(cwd)
inputfile = open(r"C:\Users\Vonag\Downloads\Project\Project\datapoints\in1.txt")
delta = 30
alphaCost = {}

def inputStringGenerator() -> list:           #reads from txt file and creates the 2 strings of input for sequence alignment
    count = 0
    inputString1 = ''
    inputString2 = ''
    with open(os.path.join(cwd, 'SampleTestCases', 'input1.txt'), 'r', encoding='UTF-8') as file:
        while line := file.readline().rstrip():
            if(line.isalpha() and not inputString1): 
                inputString1 = line
                count += 1
                continue
            elif(line.isalpha() and not inputString2):
                inputString2 = line
                count += 1
                continue

            line = int(line) + 1
            if(count == 1):
                inputString1 = inputString1[:line] + inputString1 + inputString1[line:]
            else:
                inputString2 = inputString2[:line] + inputString2 + inputString2[line:]

    return [inputString1, inputString2]


def basicSegmentNeedleman(a, b, simMatrix):
    #initialize dynamic programming array
    w, h = len(a) + 1, len(b) + 1
    #nw_array = [[0 for x in range(w)] for x in range(h)] 
    nw_array = np.zeros((w,h))
    #run the needleman wunsch algorithm to create the opt array where the opt solution is
    #initialize our array with the base cases
    for i in range(w):
        nw_array[i][0] = i*30
        
    for j in range(h):
        nw_array[0][j] = j*30

#fill out the rest of the array determinant on how the sequences match up
    for i in range(1,w):
        for j in range(1,h):
            k = a[i-1]
            l = b[j-1]
            alpha = simMatrix[k][l]
            #Match ← F(i−1, j−1) + S(Ai, Bj) where S is the similarity between the letters diagonal
            match = nw_array[i-1][j-1] + alpha
            #Delete ← F(i−1, j) + d top
            delete = nw_array[i-1][j] + 30
            #Insert ← F(i, j−1) + d left
            insert = nw_array[i][j-1] + 30
            #F(i,j) ← max(Match, Insert, Delete)
            nw_array[i][j] = min(match,delete,insert)

#CREATION OF OPT MATRIX IS CORRECT, or is it?????

#the question is do we want to create a whole other function here to create the new strings that represent
# the matches? Or should we just wrap it all up into one big function
    #i think its best to just keep it all in the same function since we can just re-use our given string lengths
    ii = len(b)
    jj = len(a)
    #create new strings with the alg now implemented
    returnString1 = ""
    returnString2 = ""
    while ii > 0 and jj > 0:
        watchValue = nw_array[ii][jj]
        match = nw_array[jj-1][ii-1] + simMatrix[a[jj-1]][b[ii-1]]
        insert = nw_array[jj][ii-1] + 30
        
        if nw_array[jj][ii] == match:
            returnString1 += a[jj-1]
            returnString2 += b[ii-1]
            ii -= 1
            jj -= 1
            
        elif nw_array[jj][ii] == insert:
            returnString1 += "_"
            returnString2 += b[ii-1]
            ii -= 1
            
        else:
            returnString1 += a[jj-1]
            returnString2 += "_"
            jj -= 1

    while ii > 0:
        returnString1 += "_"
        returnString2 += b[ii-1]
        ii -= 1
        
    while jj > 0:
        returnString1 += a[jj-1]
        returnString2 += "_"
        jj -= 1



    returnString1 = returnString1[::-1]
    returnString2 = returnString2[::-1]
    #returnString1 = "".join(reversed(returnString1))
    #returnString2 = "".join(reversed(returnString2))
    print(returnString1)
    print(returnString2)
    
    return [returnString1, returnString2, nw_array]
def divideAndConquer(a,b,simMatrix):
    
    
    return 0

def alphaMatrix():
    #set the alpha costs
    alphaCost['A'] = {"A" : 0, "C" : 110, "G" : 48, "T": 94}
    alphaCost['C'] = {"A" : 110, "C" : 0, "G" : 118, "T": 48}
    alphaCost['G'] = {"A" : 48, "C" : 118, "G" : 0, "T": 110}
    alphaCost['T'] = {"A" : 94, "C" : 48, "G" : 110, "T": 0}



if __name__ == "__main__":             #main driver
    inputs = inputStringGenerator()
    inputString1 = inputs[0]
    inputString2 = inputs[1]
    alphaMatrix()
    alpha = alphaCost
    #simMatrix = similarityMatrix(inputString1,inputString1)
    dp_array = basicSegmentNeedleman(inputString1,inputString2, alpha)
    print(inputs[0] + "\n" + inputs[1])
    #print(alpha)
    #print(dp_array)