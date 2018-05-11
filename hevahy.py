#!/usr/bin/env python3

# To import glob library, read all the pathnames matching with a specified pattern
import glob
# To import the pandas library
import pandas as pd
# To import the module 'mathplotlib.pyplot(plotting framework)'
import matplotlib.pyplot as plt
# To import numpy library as fundamental package to compute function
import numpy as np
# To import system specific parameters and functions
import sys

#Make function to count the kmers
#Define count kmers as c_kmers based on sequence and specific k
#For every i which is in range of formula (sequence length - the k value + 1),
#Define what is kmers and formula then if that kmers are not in counts then give it 0
#Otherwise, if kmers are in script, add 1 incrementally depend the number of Kmer
#Then return value from the function and pass back the counts expression

def c_kmers(seq, k):
    counts = {}
    for i in range(len(seq) - k + 1):
        kmers = seq[i:i+k]
        if kmers not in counts:
            counts[kmers] = 0
        else:
            counts[kmers] += 1
    return counts

#Define function to find other DNAs except ACGT (find the wrong DNA)
#If i is not A,C,G, or T, it is wrong DNA and add it incrementally each one
#And if the condition becomes true, output this value from this function
def others_dna(seq):
    wrong_dna = {}
    for i in seq:
        if i not in 'ACGT':
            if i in wrong_dna: wrong_dna[ i ] += 1
            else: wrong_dna[ i ] = 1
    if wrong_dna != {}: return wrong_dna

#Define function to create a table(kmers, observe, possible) in .csv
#Define data field consisting 3 columns and convert into csv
#Use index=False as default
def table_kmers(filename, data, sequename):
    tableofkmers = pd.DataFrame([[i[0], i[1], i[2]] for i in data], columns=['kmers', 'observed', 'possible'])
    tableofkmers.to_csv('result/table/'+ filename + '_' + sequename+'.csv', index=False)

#Define function to create linguistic complexity of the graph of each species name
#Name plot, x, and y label then use blue circle markers in the graph
#Give the name of plot and label then save the graph into result folder and pic name
def graph_lingcomplex(filename, sequename_l, lingcomplex):
    plt.plot(sequename_l, lingcomplex, 'bo')
    plt.xlabel('Sequence_name')
    plt.ylabel('Linguistic_Complexity')
    graph_name = filename + '_lingcomplexity.png'
    plt.savefig('result/pic/'+graph_name)

#Define function to compare between observe and possible kmers
#Make two figures (fig1 and fig2) in one graph to compare observe and possible
#Define function of graph visualization (barwidth, opacity, color, label, legend, ect)
#Save the graph into result folder and pic name
def graph_kmers(filename, sequenames, observed_kmers, possible_kmers):
    figure, axis = plt.subplots()
    index = np.arange(len(sequenames))
    barwidth = 0.2
    opacity = 0.4
    fig1 = plt.bar(index + barwidth, observed_kmers, barwidth,
                     alpha=opacity,
                     color='b',
                     label='The number of Observed Kmers')
    fig2 = plt.bar(index, possible_kmers, barwidth,
                     alpha=opacity,
                     color='r',
                     label='The number of Possible Kmers')
    plt.xlabel('Sequence')
    plt.ylabel('Kmers')
    plt.xticks(index + barwidth, sequenames)
    plt.legend()
    graph_file = filename + '_obs_poss_kmers.png'
    plt.savefig('result/pic/'+graph_file)

#Define function to create graph of wrong DNA
#Define the position and values of the graph and labels for x axis
#Make the name of the graph (.png) then save it
def wrong_dna_graph(filename, data, seq_name):
    plt.clf()
    plt.bar(range(len(data)), list(data.values()), align='center')
    plt.xticks(range(len(data)), list(data.keys()))
    graph_name = 'wrong_dna_'+ filename + '_' + seq_name + '.png'
    plt.savefig('result/pic/'+graph_name)

#Primary function
#Make the _name_ is set to the module'name
#We already import the sys module, then we need to put the filename into the script
#To work with command line arguments
#The pathnames will match the name containing *.fasta file
#Change/delete > symbol in front of sequence name
#Read Append/add the argument (sequename, total kmer etc) to the end of the list
if __name__ == "__main__":
    filename=sys.argv[1]
    if filename in glob.glob('*.fasta'):
        f = open(filename,'r')
        seq = f.readlines()
        filename = filename.replace('.', '_')
        sequenames = []
        observed_kmers = []
        possible_kmers = []
        lingcomplex = []
        for line_num, line in enumerate(seq[0:len(seq)]):
            if len(line) > 1 :
                if '>' in line :
                    line = line.replace(">", "")
                    sequename = line.rstrip()
                    sequenames.append(sequename)
                else:
                    seq = line.rstrip()
                    #check per each sequence name, the sequence format
                    seque_format = others_dna(seq)
                    klist = []
                    possibles = []
                    observes = []
                    for k in range(1,len(seq)+1):
                        #check the possible kmers
                        if 4**k < len(seq):
                            possible = 4**k
                        else:
                            possible = len(seq) - k + 1
                        #get the kmers
                        counts = c_kmers(seq, k)
                        #print(counts)
                        #get the observed kmers
                        observed = len(counts)
                        klist.append(k)
                        possibles.append(possible)
                        observes.append(observed)
                    #get the total possible kmers
                    possible_total = sum(possibles)
                    possibles.append(possible_total)
                    possible_kmers.append(possible_total)
                    #get the total observed kmers
                    observed_total = sum(observes)
                    observes.append(observed_total)
                    observed_kmers.append(observed_total)
                    klist.append('Total Kmers');
                    #Combine data of klist, observes and possible
                    data = list(zip(klist, observes, possibles))
                    #Create table consisting the data and sequence name of the filename
                    det = table_kmers(filename, data, sequename)
                    #Define function the linguistic complexity then append the argument
                    linguistic_complexity = observed_total/possible_total
                    lingcomplex.append(linguistic_complexity)
        #Print all linguistic complexity
        #Make the language of execution
        all_lingcomplex = graph_lingcomplex(filename, sequenames, lingcomplex)
        if seque_format != None:
            wrong_graph = wrong_dna_graph(filename, seque_format, sequename)
            print('Wrong dna exist in', sequename, ' detail:', seque_format)
        summary_kmers = graph_kmers(filename, sequenames, observed_kmers, possible_kmers)
        print('Success')

    else:
        print('Not success')
