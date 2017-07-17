from random import randint
import numpy as np
import pandas as pd
import random
from gensim.models import word2vec
from sklearn.cluster import KMeans
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
<<<<<<< HEAD
from sklearn.metrics import roc_auc_score
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from sklearn.metrics import precision_score
from sklearn.metrics import accuracy_score
from sklearn.metrics import recall_score
from xgboost import XGBClassifier
=======
from sklearn.metrics import precision_recall_fscore_support
from sklearn.metrics import roc_auc_score
from sklearn.svm import SVC
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.feature_extraction.text import TfidfVectorizer
from collections import Counter, defaultdict
from sklearn.pipeline import Pipeline
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.metrics import precision_score
from sklearn.metrics import accuracy_score
from sklearn.metrics import recall_score
>>>>>>> origin/master

def windower(sequence, position, wing_size):
    # window size = wing_size*2 +1
    position = int(position)
    wing_size = int(wing_size)
    if (position - wing_size) < 0:
        return sequence[:wing_size + position]
    if (position + wing_size) > len(sequence):
        return sequence[position - wing_size:]
    else:
        return sequence[position - wing_size:position + wing_size]

class DataCleaner:
<<<<<<< HEAD
    def __init__(self, output, data="phosphosites.csv", delimit=",", amino_acid="K", window_size=7, training_ratio=.7, header_line=0, neg_per_seq=2):
=======
    def __init__(self, output, data="phosphosites.csv", delimit=",", amino_acid="K", sites="code",
                 modification="phosphorylation", window_size=7, pos="position", training_ratio=.7,
                 header_line=0, seq="sequence", neg_per_seq=2, lines_to_read=10000):
>>>>>>> origin/master

        data = pd.read_csv(data, header=header_line, delimiter=delimit, quoting=3, dtype=object)
        self.data = data.reindex(np.random.permutation(data.index))
        self.amino_acid = amino_acid
        self.training_ratio = training_ratio  # Float value representing % of data used for training
        self.proteins = {}
        self.neg_count = 0
        self.neg_per_seq = neg_per_seq
        self.window = int(window_size)
        self.features= []
        self.labels = []
        self.output = open(output, "a")
        sequences = self.data["sequence"]
        positive_sites = self.data["position"]
        size = len(self.data["sequence"])
        for i in range(0,size):
            #print(sequences[i][int(positive_sites[i])-1])
            try:
                self.features.append(windower(sequences[i], positive_sites[i], self.window))
                self.labels.append(1)
            except:
                print(i)

        counter = len(self.features)
        for i in range(int(counter*neg_per_seq)):
            if len(self.features) >= counter*neg_per_seq:
                break
            selector = randint(0, size)
            options = []
            try:
                for j in range(len(sequences[selector])):
                    if sequences[selector][j] == self.amino_acid:
                        options.append(j)
            except:
                pass
            if len(options) > 0:
                try:
                    random.shuffle(options)
                    for j in options:
                        t = windower(sequences[selector],j,self.window)
                        if t not in self.features:
                            self.features.append(t)
                            self.labels.append(0)
                except:
                    pass


        temp = list(zip(self.features, self.labels))
        random.shuffle(temp)
        self.features, self.labels = zip(*temp)
        print(len(self.features), len(self.labels))
        for i in range(len(self.features)):
            t = str(self.features[i])+","+str(self.labels[i])+"\n"
            self.output.write(t)

class Classy:
<<<<<<< HEAD
    def __init__(self, benchmark="benchmark.csv", data="clean_serine.csv", delimit=",", amino_acid="Y", training_ratio=.7, header_line=0):
=======
    def __init__(self, data="clean_serine.csv", delimit=",", amino_acid="Y", training_ratio=.7, header_line=0):
>>>>>>> origin/master
        self.data = open(data, "r")
        self.amino_acid = amino_acid
        self.training_ratio = training_ratio  # Float value representing % of data used for training
        self.features= []
        self.labels = []
        i = 0
        for line in self.data:
            try:
                x, y = line.split(",")
                y = int(y.strip("\n"))
                t = []
                for j in x:
                    t.append(j)
                self.features.append(t)
                self.labels.append(y)
            except:
                print("Bad data at line"+str(i))
            i = i + 1
        temp = list(zip(self.features, self.labels))
        random.shuffle(temp)

        self.features, self.labels = zip(*temp)
<<<<<<< HEAD
        self.features = list(self.features)
        self.labels = list(self.labels)
        benchmark = pd.read_csv(benchmark,header=header_line, delimiter=delimit, quoting=3, dtype=object)
        s =  benchmark["fasta_seq"]
        a = benchmark["modification_residue"]
        p = benchmark["modification_region_location"]
        for i in range(len(a)):
            if a[i] == self.amino_acid:
                self.features.append(s[i])
                self.labels.append(1)


        self.num_features = 200  # Word vector dimensionality
=======

        self.num_features = 300  # Word vector dimensionality
>>>>>>> origin/master
        self.min_word_count = 1  # Minimum word count
        self.num_workers = 4  # Number of threads to run in parallel
        self.context = 5  # Context window size
        self.downsampling = 5e-1  # Downsample setting for frequent words
        self.model = word2vec.Word2Vec(self.features ,workers=self.num_workers, size=self.num_features, min_count=self.min_word_count,window=self.context, sample=self.downsampling)



    def kluster(self):
        word_vectors = self.model.wv.syn0
<<<<<<< HEAD
        num_clusters = 15
=======
        num_clusters = 15 # og is 4
>>>>>>> origin/master
        print(num_clusters)
        kmeans_clustering = KMeans(n_clusters=num_clusters)
        idx = kmeans_clustering.fit_predict(word_vectors)
        word_centroid_map = dict(zip(self.model.wv.index2word, idx))
        for cluster in range(0, 10):
            print("Cluster" +str(cluster))
            words = []
            val = list(word_centroid_map.values())
            key = list(word_centroid_map.keys())
            for i in range(len(val)):
                if val[i] == cluster:
                    words.append(key[i])
            print(words)
            train_centroids = np.zeros((len(self.features), num_clusters),dtype="float32")
            counter = 0
            for sequence in self.features:
                train_centroids[counter] = bag_of_centroids(sequence,  word_centroid_map)
                counter += 1
<<<<<<< HEAD
        self.features = list(self.features)
        self.labels = list(self.labels)
        print(len(train_centroids),len(self.features))

        t1 = train_centroids

        for i in range(len(self.features)):
            try:
                seq = "".join(self.features[i])
                t = ProteinAnalysis(seq).gravy()
                np.append(t1[i], t)
            except:
                pass

        train_centroids = t1

        s = len(train_centroids)
        test_centroids = train_centroids[int(.9*s):]
        self.test_labels = self.labels[int(.9*s):]
        train_centroids = train_centroids[:int(.9*s)]
        self.labels = self.labels[:int(.9*s)]
        X_train, X_test, y_train, y_test = train_test_split(train_centroids, self.labels, test_size = 0.1, random_state = 42)
        np.append(X_test, test_centroids, axis=0)
        np.append(y_test, self.test_labels, axis=0)
        forest = RandomForestClassifier(n_estimators=100)
        #forest = XGBClassifier()
        forest.fit(X_train, y_train)
        result = forest.predict(X_test)

        custom_results(y_test, result)
=======
        X_train, X_test, y_train, y_test = train_test_split(train_centroids, self.labels, test_size = 0.33, random_state = 42)
        forest = RandomForestClassifier(n_estimators=100)

        forest.fit(X_train, y_train)
        result = forest.predict(X_test)

        print(precision_score(y_test, result))
        print(recall_score(y_test, result))
        print(accuracy_score(y_test, result))
>>>>>>> origin/master
        print(roc_auc_score(y_test, result))


def bag_of_centroids(wordlist, word_centroid_map):
    num_centroids = max(word_centroid_map.values()) + 1
    bag_of_centroids = np.zeros(num_centroids, dtype="float32")
    for word in wordlist:
        if word in word_centroid_map:
            index = word_centroid_map[word]
            bag_of_centroids[index] += 1
    return bag_of_centroids

<<<<<<< HEAD
def custom_results( true_y, y):
    tp = 0
    tn = 0
    fn = 0
    fp = 0
    for i in range(len(y)):
        if y[i] == true_y[i]:
            if true_y[i] == 1:
                tp = tp + 1
            else:
                tn = tn + 1
        else:
            if y[i] == 1:
                fp = fp + 1
            else:
                fn = fn + 1
    recall = tp / (tp + fn)
    precision = tp / (tp+fp)
    accuracy = (tp + tn) / (tp + fp + tn +fn)
    print(precision, accuracy, recall)
#y= DataCleaner(amino_acid="K", data="k_site.csv", output="clean_k.csv")

x= Classy(data="clean_s.csv", amino_acid="S")
=======

y= DataCleaner(amino_acid="K", data="k_site.csv", output="clean_k.csv")
x= Classy(data="clean_k.csv")
>>>>>>> origin/master
x.kluster()
