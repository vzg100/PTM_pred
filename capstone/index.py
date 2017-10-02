import numpy as np
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import precision_recall_fscore_support
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import random
from random import randint
import numpy as np
from matplotlib import pyplot as plt

def windower(seq, pos, w):
    if (pos - w) < 0:
        return seq[:w+pos]
    if (pos + w) > len(seq):
        return seq[pos-w:]
    else:
        return seq[pos-w:pos+w]

class Classy:
    
    def __init__(self, data="phosphosites.txt", delimit=",", amino_acid="Y", sites="code",
                 modification="phosphorylation", window_size=7, pos="position", training_ratio=.7,
                 header_line=0, seq="sequence", neg_per_seq = 10, lines_to_read=50000, forest_size=110):
        
        data = pd.read_csv(data, header=header_line, delimiter=delimit, quoting=3)
        self.data = data.reindex(np.random.permutation(data.index))
        self.amino_acid = amino_acid
        self.training_ratio = training_ratio
        self.protiens = {}
        self.features = []
        self.labels = []
        for i in range(lines_to_read):
            if("X" not in data[seq][i]) and (data[sites][i] == amino_acid) and (data[seq][i] not in self.protiens.keys()):
                self.protiens[data[seq][i]]=[data[pos][i]]
            elif("X" not in data[seq][i]) and (data[sites][i] == amino_acid) and (data[pos][i] not in self.protiens[data[seq][i]]):
                self.protiens[data[seq][i]].append(data[pos][i])
            else:
                pass
        #neg_per_seq = len(self.protiens.keys())
        for i in self.protiens.keys():
            neg_pos = []
            neg_pos_used = []
            
            positions = self.protiens[i]
            for p in positions:
                window = ProteinAnalysis(windower(i, p, window_size))
                self.features.append([window.gravy(), window.aromaticity(), window.isoelectric_point()])
                self.labels.append(1)
            for aa in range(len(i)):
                if (i[aa] == amino_acid) and (aa not in positions):
                    neg_pos.append(aa)
            counter = 0
            while(counter < neg_per_seq) and (len(neg_pos_used) != len(neg_pos)):
                qq = randint(0, len(neg_pos))
                if(qq not in neg_pos_used):
                    window = ProteinAnalysis(windower(i, qq, window_size))
                    self.features.append([window.gravy(), window.aromaticity(), window.isoelectric_point()])
                    self.labels.append(0)   
                    
                    neg_pos_used.append(qq)
                    counter += 1

        temp = list(zip(self.features, self.labels))
        random.shuffle(temp)
        self.features, self.labels = zip(*temp) 
    def calculate(self):
        training_slice = int(self.training_ratio*len(self.labels))
        self.forest = RandomForestClassifier(n_estimators = forest_size, verbose = 0, n_jobs = 4)
        self.forest = self.forest.fit(self.features[0:training_slice], self.labels[0:training_slice])
        self.results = self.forest.predict(self.features[training_slice:])
        self.rating = precision_recall_fscore_support(self.labels[training_slice:], self.results, average="macro")
        
        
    def scores(self):
        tp = 0
        tn = 0
        fp = 0
        fn = 0
        for i in range(len(self.results)):
            if self.results[i] == 1 and self.labels[i+training_slice] == 1:
                tp+=1
            elif self.results[i] == 0 and self.labels[i+training_slice] == 0:
                tn +=1 
            elif self.results[i] == 1 and self.labels[i+training_slice] == 0:
                fp +=1
            else:
                fn+=1
        self.prec = tp/(tp+fp) # sensitivity 
        self.recall = tp/(tp+fn)
        self.acc = (tp+tn)/(tp+fp+tn+fn)
        print(self.prec, self.recall, self.acc)
        #count positives vs negatives and where they pop up, look for feedback places online 
        #take run of the mill problem and compare to accruacy 
        #check cbs paper 
    def save_to_csv(self, filname):
        f = open(filename, "w")
        f.write(self.amino_acid+","+self.prec+","+self.recall+","+self.acc+"\n")
        