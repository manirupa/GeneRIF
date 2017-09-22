#!/usr/bin/python

import time
import operator
import sys, re, os
from math import floor
from random import randint

# gensim modules
from gensim import utils
from gensim.models.doc2vec import LabeledSentence
#from gensim.models import Doc2Vec
from gensim.models import *

# numpy
import numpy

# random
from random import shuffle

# classifier
from sklearn.linear_model import LogisticRegression


ENT_RE = re.compile(r'[,;.(){}]')

def remove_entities(text):
    text = text.lower()
    return ENT_RE.sub('', text)

def introspect(desc, list):
    print desc
    for i in range(len(list)):
        print 'Record %s: %s\n' % (i+1, list[i])
    print 'Length ', len(list)

def write_to_file(str_to_print,filename):
    f = open(filename, 'ab')
    f.write(str_to_print)
    f.close()
    return

def check():
    if (len(sys.argv)<4):
        print "usage: %s <model_params>: gene_id, size, window [optional: <gramsize>]"  % sys.argv[0]
        exit()
    return

def write_to_file(str_to_print,filename):
    f = open(filename, 'ab')
    f.write(str_to_print)
    f.close()
    
def format_vec(id, vec):
    ftdvec = id[0]
    for v in vec:
        ftdvec = '%s,%s' % (ftdvec,v)
    return ftdvec
    
def write_model(model,modelname,suffix,docs,tags):
    filename = '%s_%s' % (modelname,suffix)
    n = len(docs)
    for i in range(n):
        docvec = model.docvecs[i]
        str_to_print = '%s\n' % format_vec(tags[i],docvec)
        #str_to_print = '%s\t%s\n' % (docs[i][1],docvec)
        #print format_vec(docs[i][1],docvec)
        write_to_file(str_to_print,filename)
    return

def do_this():
    check()
    gene_id = sys.argv[1]
    size = sys.argv[2]
    window = sys.argv[3]
    gramsize = 1
    if(len(sys.argv)>4):
        gramsize = sys.argv[4]
    if(int(gramsize) > 1):
        #e.g.model_s100_w4_g26_generif10-50K.3-gram-cluster.txt
        model_cfile = 'model_s%s_w%s_g%s_generif10-50K.%s-gram-cluster.txt' % (size,window,gene_id,gramsize)
    else:
        model_cfile = 'model_s%s_w%s_g%s_generif10-50K-cluster.txt' % (size,window,gene_id)
    genefile = 'GeneDatafilesTM/Gene%s.txt' % (gene_id)
    ncbifile = 'GeneDatafilesTM/NCBI%s.txt' % (gene_id)
    print "Cluster file: ", model_cfile
    print "Gene file: ", genefile
    print "NCBI file: ", ncbifile
    
    clines = open(model_cfile).readlines()
    glines = open(genefile).readlines()
    nlines = open(ncbifile).readlines()
    
    generifs = []
    for i in range(len(glines)):
        generifs.append(glines[i].split('\t')[4].strip())
    print generifs
    
    ncbi = []
    for i in range(len(nlines)):
        ncbi.append(nlines[i].strip())
    print ncbi

    #Process clusters, from cluster file,gene file 
    clusters = {}

    #Make clusters
    for i in range(len(clines)):
        c_id = int(clines[i])
        try:
            (clusters[c_id]).append(generifs[i])
        except:
            clusters[c_id] = []
            (clusters[c_id]).append(generifs[i])

    #print clusters   
    
    #for i in clusters.keys():
    #    print 'cluster %s: \n' % i, clusters[i]
    
    orig_stdout = sys.stdout
    ranking_file = model_cfile.replace('cluster','ranking')
    f = file(ranking_file, 'w')
    sys.stdout = f
    tic = time.time()
    #Gen Rankings
    rankings = []
    for i in clusters.keys():
        length = len(clusters[i])
        index = int(floor(length/2))#randint(0,length-1)
        print "cluster: ", i, "length: ", length, "index: ",  index
        rankings.append(clusters[i][index])
    
    print "Rankings:", rankings
    print "NCBI:", ' '.join([x for x in ncbi])
    print "Rankings:", ' '.join([x for x in rankings])

    toc = time.time()
    te = toc - tic
    print "Time elapsed for training model: ", te 
    sys.stdout = orig_stdout
    f.close()
    
    return

do_this()