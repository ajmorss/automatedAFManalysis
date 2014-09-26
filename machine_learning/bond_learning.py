# -*- coding: utf-8 -*-

import csv
import numpy as np
import random
from sklearn import svm, metrics
from sklearn import cross_validation
import cPickle

class FeatureInfo(object):
    #Object that calculates and stores the means and stds of the various
    #features
    def __init__(self, features):
        self.featuremeans = []
        self.featurestds = []
        for f in xrange(0, features.shape[1]):
            self.featuremeans.append(np.mean(features[:,f]))
            self.featurestds.append(np.std(features[:,f]))

def shuffletrain(samples, targets):
    #Randomizes the order of the training data
    c=[]
    for x in xrange(0,len(targets)):
        temp=[]
        for y in xrange(0,len(samples[0])):
            try:
                temp.append(samples[x][y])
            except IndexError:
                pass
        temp.append(targets[x])
        c.append(temp)
    random.shuffle(c)
    sampout=[]
    targout=[]
    for x in c:
        targout.append(x[len(x)-1])
        sampout.append(x[0:len(x)-1])
    return sampout,targout

def cusscore(ground_truth,predictions):
    #Custom sklearn scoring function
    mat = metrics.confusion_matrix(ground_truth, predictions)
    print mat
    return (float(mat[1,1])/(float(mat[0,1])+float(mat[1,1])+float(mat[2,1]))**2
            *float(mat[1,1])/(float(mat[1,0])+float(mat[1,1])+float(mat[2,0])))

def builddat():
    #Loads features and classifications from file, randomizes their order,
    #finds mean and std then adjusts them so their mean is 0 and their STD is 1
    features = []
    targets = []
    with open('bondsamples.txt','r') as datin:
        datin = csv.reader(datin, delimiter=',')
        for rows in datin:
            try:
                for y in rows:
                    if np.isnan(float(y)):
                        raise ValueError
                    if np.isinf(float(y)):
                        raise ValueError
                temp = [float(y) for y in rows]
                if len(temp)==41:
                    features.append(temp[:-2:])
                    targets.append(temp[-2])
            except ValueError:
                pass

    features = np.array(features)
    targets = np.array(targets)
    (features, targets) = shuffletrain(features, targets)
    features = np.array(features)
    targets = np.array(targets)
    featureinfo = FeatureInfo(features)
    cPickle.dump( featureinfo, open( "feature_info_for_bond_classifier2.p", "wb" ) )
    for f in xrange(0, features.shape[1]):
        if np.std(features[:,f])==0:
            features[:,f] = (features[:,f]-np.mean(features[:,f]))/1
        else:
            features[:,f] = (features[:,f]-np.mean(features[:,f]))/np.std(features[:,f])
    return features, targets

#Entry point
if __name__ == '__main__':

    #Make the scorer, load the data
    my_custom_scorer = metrics.make_scorer(cusscore, greater_is_better=True)
    features, targets_dat = builddat()
    targets_dat = np.ravel(targets_dat)

    #RBF SVM Optimization
    best = [0,0,0]
    #Insert a list of parameters to test over here, carr is a list of C values
    #and garr is a list of gamma values for the sklearn SVC
    carr = [20]
    garr = [0.6]
    print targets_dat.size
    #For all values of carr, garr train a classifer with crossvalidation
    for c in carr:
        for g in garr:
            classifier = svm.SVC(C = c, gamma = g, class_weight='auto')
            scores = cross_validation.cross_val_score(classifier, features,
                                                      targets_dat, cv=5,
                                                      scoring=my_custom_scorer)
            #Keep the parameters of the classifier that gives the best score
            if best[0]<np.mean(scores):
                best[0]=np.mean(scores)
                best[1]=c
                best[2]=g

    #Final RBF SVM Training
    classifier = svm.SVC(C=best[1],gamma=best[2],class_weight='auto')
    classifier.fit(features, targets_dat)
    predictions = classifier.predict(features)
    print metrics.confusion_matrix(targets_dat, predictions)
    cPickle.dump( classifier, open( "best_tether_classifier2.p", "wb" ) )
