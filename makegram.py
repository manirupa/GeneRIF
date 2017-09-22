#!/usr/bin/python

"""
Author: Manirupa Das
This script makes n-gram data files given original data files and input n
"""

import time
import operator
import sys, re, os
from textblob import *
from random import randint

# gensim modules
from gensim import utils
from gensim.models.doc2vec import LabeledSentence
#from gensim.models import Doc2Vec
from gensim.models import *

from nltk.corpus import stopwords
stop = stopwords.words('english')

# numpy
import numpy

# random
from random import shuffle

# classifier
from sklearn.linear_model import LogisticRegression


ENT_RE = re.compile(r'[,;.(){}-]')

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
        print "usage: %s <input_file> <field_num> <N for n-gram>"  % sys.argv[0]
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

def make_gram(sentence, gramsize):
    s_grams = []
    text = remove_entities(sentence)
    temp = text.split()
    n = len(temp)
    for i in range(n-(gramsize-1)):
        gram = '-'.join([x for x in temp[i:i+(gramsize)]])
        if (i%1000 == 0):
            print i, gram
        s_grams.append(gram)
    final_grams = ' '.join([s for s in s_grams])
    print final_grams
    return final_grams
    

def do_this():
    check()
    input_file = sys.argv[1]
    fieldnum = int(sys.argv[2])
    gramsize = int(sys.argv[3])

    print "Gene file: ", input_file
    glines = open(input_file).readlines()
    
    generifs = []
    for i in range(len(glines)):
        temp = glines[i].split('\t')
        text = temp[fieldnum].strip()
        rifgrams = make_gram(text,gramsize)
        temp[fieldnum] = rifgrams
        generifs.append(temp)
        outfile = '%s.%s.txt' % (input_file.replace('.txt',''), '%s-gram' % gramsize)
        str_to_print = '%s\n' % '\t'.join([t for t in temp]) 
        write_to_file(str_to_print, outfile)
    print generifs
    
    return

do_this()