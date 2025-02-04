"""
Protein Sequencing Project
Name:
Roll Number:
"""

from os import X_OK
import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    f=open(filename)
    x=f.read().splitlines()
    y="".join(x)
    return y


'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    x=[]
    y=[]
    for i in range(startIndex, len(dna), 3):
        x.append(dna[i:i+3])
        if dna[i:i+3]=='TAG' or dna[i:i+3]=='TAA' or dna[i:i+3]=='TGA':
            break
    for string in x:
        string=string.replace("T","U")
        y.append(string)
    return y


'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):
    import json
    a={}
    f=open(filename)
    b=json.load(f)
    for x,y in b.items():
        for i in y:
            c=i.replace("T", "U")
            a[c]=x
    return a


'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):
    x=[]
    if codons[0]=='AUG':
        x.append("Start")
        for i in range(1,len(codons)):
            if codons[i] in codonD.keys():
                x.append(codonD[codons[i]])
        return x


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    x=readFile(dnaFilename)
    y=makeCodonDictionary(codonFilename)
    i=0
    j=0
    k=[]
    while i<len(x):
        if x[i:i+3]=="ATG":
            m=dnaToRna(x,i)
            n=generateProtein(m,y)
            k.append(n)
            i=i+3*len(m)
        else:
            i+=1
            j+=1
    return k


def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    m=[]
    for i in proteinList1:
        for j in proteinList2:
            if i==j and i not in m:
                m.append(i)
    return m


'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    a=[]
    for i in proteinList:
        for j in i:
            if i not in a:
                a.append(j)
    return a


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    x={}
    for i in aaList:
        if i not in x:
            x[i]=1
        else:
            x[i]+=1
    return x


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    x=combineProteins(proteinList1)
    y=combineProteins(proteinList2)
    p=aminoAcidDictionary(x)
    q=aminoAcidDictionary(y)
    t=[]
    diff=[]         
    m={}
    n={}    
    for i in p:
        m[i] = p[i]/len(x)
        if i not in t and i !="Start" and i !="Stop":
            t.append(i)
        # print(m[i])
    for j in q:
        n[j] = q[j]/len(y)
        if j not in t and j !="Start" and j!="Stop":
            t.append(j)
    for k in t:
        f1=0
        f2=0
        if k in m:
            f1= m[k]
        if k in n:
            f2= n[k]
        difference = f2-f1
        if difference < -cutoff or difference > cutoff:
            z=[k, f1, f2]
            diff.append(z)
    return diff


'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    print("The following proteins occurred in both DNA Sequences:")
    for i in commonalities:
        x=""
        l=i[1:(len(i)-1)]
        c=0
        for j in l:
            x+=j
            c+=1
            if c!=len(l):
                x+="-"
        if len(x)!=0:
            print(x)
    print("The following amino acids occurred at very different rates in the two DNA sequences:" )
    for i in differences:
        m=i[0]
        f1=round(i[1]*100,2)
        f2=round(i[2]*100,2)
        print(str(m)+" "+str(f1)+" % in Seq1"+","+str(f2)+"% in Seq2")
    return


def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    l=[]
    x=combineProteins(proteinList1)
    y=combineProteins(proteinList2)
    dict1=aminoAcidDictionary(x)
    dict2=aminoAcidDictionary(y)
    for i in dict1:
        if i not in l:
            l.append(i)
    for j in dict2:
        if j not in l:
            l.append(j)
    l.sort()
    return l


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    m=combineProteins(proteinList)
    n=aminoAcidDictionary(m)
    x=[]
    for i in labels:
        if i in n:
            x.append(n[i]/len(m))
        else:
            x.append(0)
    return x


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    import numpy as np
    w=0.4
    x=np.arange(len(xLabels))
    y=np.arange(-w,len(xLabels)-1,1)
    plt.bar(y,freqList1,width=-w,label=label1,edgecolor=edgeList)
    plt.bar(x,freqList2,width=w,label=label2,edgecolor=edgeList)
    plt.xticks(ticks=list(range(len(xLabels))),labels=xLabels,rotation="vertical")
    plt.legend()
    plt.title("Comparision of Frequencies")
    plt.show()
    return


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    x=[]
    y=[]
    for i in range(len(biggestDiffs)):
        y.append(biggestDiffs[i][0])
    for i in range(len(labels)):
        if labels[i] in y:
            x.append("black")
        else:
            x.append("white")
    return x


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():
    human_proteins=synthesizeProteins("data/human_p53.txt","data/codon_table.json")
    elephant_proteins=synthesizeProteins("data/elephant_p53.txt","data/codon_table.json")
    pro=commonProteins(human_proteins,elephant_proteins)
    diff=findAminoAcidDifferences(human_proteins,elephant_proteins,0.005)
    displayTextResults(pro,diff)
    labels=makeAminoAcidLabels(human_proteins,elephant_proteins)
    f1=setupChartData(labels,human_proteins)
    f2=setupChartData(labels,elephant_proteins)
    edges=makeEdgeList(labels,diff)
    createChart(labels, f1, "Human", f2, "Elephant", edgeList=edges)
    return


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    # print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    # test.week1Tests()
    # print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    # runWeek1()
    # # test.testReadFile()
    # # test.testDnaToRna()
    # # test.testMakeCodonDictionary()
    # # test.testGenerateProtein()
    # test.testSynthesizeProteins()

    # ##Uncomment these for Week 2 ##
    
    # print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    # test.week2Tests()
    # print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    # runWeek2()
    
    # # test.testCommonProteins()
    # # test.testCombineProteins()
    # # test.testAminoAcidDictionary()
    # # test.testFindAminoAcidDifferences()

    # ## Uncomment these for Week 3 ##
    
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()

    # test.testMakeAminoAcidLabels()
    # test.testSetupChartData()