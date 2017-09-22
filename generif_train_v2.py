#!/usr/bin/python

import time
import operator
import sys, re, os
from textblob import *

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
        print "usage: %s <input_file> <size> <window>" % sys.argv[0]
        exit()
    return

def process_data(lines):
    data = []
    documents = []
    for i in range(len(lines)):
        temp = lines[i].split('\t')
        try:
            taxid = temp[0].strip()
        except:
            taxid = 't999999'
        try:
            geneid = temp[1].strip()
        except:
            geneid = 'g999999'
        try: 
            pubmedid = temp[2].strip()
        except:
            pubmedid = 'p999999'
        try:
            text = temp[4].strip()
        except:
            text = 'Text999999'
        text = remove_entities(text)
        label = 'SENT_%s-%s-%s' % (pubmedid,taxid, geneid)
        thewords = text.split()
        data.append([thewords,[label]])
        lbldsent = LabeledSentence(words=thewords, tags=[label,pubmedid,taxid,geneid])
        documents.append(lbldsent)
    #introspect('generif data', data[200:210])
    return [data,documents]

def get_documents(data):
    _words = []
    _labels = []
    for d in data:
        _words.append(d[0])
        _labels.append(d[1])
    documents = LabeledSentence(_words,_labels) #create in doc2vec format #prints but wrong
    return documents

def write_to_file(str_to_print,filename):
    f = open(filename, 'ab')
    f.write(str_to_print)
    f.close()
    
def format_vec(id, vec):
    ftdvec = id[0]
    for v in vec:
        ftdvec = '%s,%s' % (ftdvec,v)
    return ftdvec
    
def write_model(model,modelname,suffix,data):
    filename = '%s_%s' % (modelname,suffix)
    docs = data[0]
    n = len(docs)
    for i in range(n):
        docvec = model.docvecs[i]
        str_to_print = '%s\n' % format_vec(docs[i][1],docvec)
        #str_to_print = '%s\t%s\n' % (docs[i][1],docvec)
        #print format_vec(docs[i][1],docvec)
        write_to_file(str_to_print,filename)
    return

def do_this():
    check()
    infile = sys.argv[1]
    size = sys.argv[2]
    window = sys.argv[3]
    lines = open(infile).readlines()
    data = process_data(lines)
    documents = data[1]

    #test structure
    #introspect('generif docs', documents[0:10])
    #introspect('data', data[0:10])
    orig_stdout = sys.stdout
    suffix = os.path.basename(infile)
    f = file('model_s%s_w%s_%s.log.txt' % (size,window,suffix), 'w')
    sys.stdout = f
    print "Model params: size - %s, window - %s, corpus - %s" % (size,window,suffix)
    print "length data: ", len(data), " length docs:", len(documents)
    print "length data[0]: ", len(data[0]), " length data[1]:", len(data[1])

    #train different models
    #1. size 100
    tic = time.time()
    model = doc2vec.Doc2Vec(documents, min_count = 1, window = window, size = size, workers=20, dm_tag_count = 4, dbow_words = 1) #dm_concat = 1,
    toc = time.time()
    te = toc - tic
    print "Time elapsed for training model: ", te 
    n = 999
    docvec = model.docvecs[n] 
    #print 'docvec  %s: ' % n , docvec, len(docvec)
    print "Model sim:", model.docvecs.most_similar("SENT_22012226-9606-10"), len(model.docvecs.most_similar("SENT_22012226-9606-10"))
    write_model(model,'model_s%s_w%s' % (size,window),suffix,data)
    
    sys.stdout = orig_stdout
    f.close()
    return

do_this()