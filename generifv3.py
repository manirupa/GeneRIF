#!/usr/local/bin/python

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

def check():
    if (len(sys.argv)<5):
        print "usage: %s <input_file> <size> <window> <gene_id_list_file>" % sys.argv[0]
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
        write_to_file(str_to_print,filename)
    return

def process_data(lines, geneids):
    data = []
    documents = []
    #also makes sets by geneid
    sets = {}
    tags = {}

    for i in range(len(lines)):
        temp = lines[i].split('\t')
        #print temp
        try:
            taxid = temp[0].strip()
        except:
            taxid = 't999999'
            print 'taxid', i, lines[i]
        try:
            geneid = temp[1].strip()
        except:
            geneid = '999999'
            print 'geneid', i, lines[i]
        try: 
            pubmedid = temp[2].strip()
        except:
            pubmedid = 'p999999'
            print 'pubmedid', i, lines[i]
        try:
            text = temp[4].strip()
        except:
            text = 'Text999999'
            print 'text',i, lines[i]
        text = remove_entities(text)
        label = 'SENT_%s-%s-%s' % (pubmedid,taxid, geneid)
        thewords = text.split()
        data.append([thewords,[label]])
        lbldsent = LabeledSentence(words=thewords, tags=[label,pubmedid,taxid,geneid])
        geneid = int(geneid)
        try:
            (sets[geneid]).append(lbldsent)            
        except:
            #print "yoohoo"
            sets[geneid] = []
            (sets[geneid]).append(lbldsent)            
        try:
            (tags[geneid]).append([label,pubmedid,taxid,geneid])
        except:
            #print "yoohoo"
            tags[geneid] = []
            (tags[geneid]).append([label,pubmedid,taxid,geneid])
        documents.append(lbldsent)
    #introspect('generif data', data[200:210])
    return [data,documents,sets,tags]


def do_this():
    check()
    infile = sys.argv[1]
    size = sys.argv[2]
    window = sys.argv[3]
    geneids = eval(open(sys.argv[4]).read())
    #print geneids
    lines = open(infile).readlines()
    data = process_data(lines, geneids)
    documents = data[1]
    sets = data[2]
    tags = data[3]

    print "num sets: ", len(sets.keys())
    #print "sets:", sets
    #print sets[9], len(sets[9])
    doc = sets[9][1]
    tagsfor9 = tags[9]
    print "doc:", doc
    print "tags for 9:", tagsfor9
   
    #test structure
    #print 'generif docs: ', documents[0:5] 
    orig_stdout = sys.stdout
    suffix = os.path.basename(infile)
    model_folder = 'model_s%s_w%s_%s' % (size,window,suffix.replace('.txt',''))
    os.system('mkdir %s' % model_folder)
    f = file('%s/model_s%s_w%s_%s.log.txt' % (model_folder,size,window,suffix), 'w')
    sys.stdout = f
    print "Model params: size - %s, window - %s, corpus - %s" % (size,window,suffix)
    print "length data: ", len(data), " length docs:", len(documents)
    print "length data[0]: ", len(data[0]), " length data[1]:", len(data[1])

    #train different models
    #1. size 100
    tic = time.time()
    for g in geneids:
        print g
        docset = sets[g]
        tagset = tags[g]
        print 'gene: ', g, len(tagset)
        model = doc2vec.Doc2Vec(docset, min_count = 1, window = int(window), size = int(size), workers=20, \
                                dm_concat = 1,  dm_tag_count = 4 , dbow_words = 1) 
        #n = 1
        #docvec = model.docvecs[n] 
        #print 'docvec  %s: ' % n , docvec, len(docvec)
        #print "Model sim:", model.docvecs.most_similar("SENT_22012226-9606-10"), len(model.docvecs.most_similar("SENT_22012226-9606-10"))
        modelfile = '%s/%s' % (model_folder,'model_s%s_w%s_g%s' % (size,window,g))
        write_model(model,modelfile,suffix,docset,tagset)
    
    toc = time.time()
    te = toc - tic
    print "Time elapsed for training model: ", te 
    sys.stdout = orig_stdout
    f.close()
    
    return

do_this()
