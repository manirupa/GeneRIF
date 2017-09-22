#!/usr/local/bin/python

"""
@author: Manirupa Das
This script normalizes data, and given experimental settings, performs 
doc2vec training on sets of GenRIFs on a gene-by-gene basis (within-gene)
"""

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
    text = text.lower() #already lower-casing
    return ENT_RE.sub('', text)

def introspect(desc, list):
    print desc
    for i in range(len(list)):
        print 'Record %s: %s\n' % (i+1, list[i])
    print 'Length ', len(list)

def check():
    if (len(sys.argv)<4):
        print "usage: %s <input_file> <gene_id_`list format`_file> <experimental_setting_`list format`_file>" % sys.argv[0]
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
    uniqtext = []
    
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
        record = '%s %s' % (geneid,text)
        #check for dupes, if not dupe then add
        if (record not in uniqtext):
        	uniqtext.append(record)
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
        else:
            pass
    #introspect('generif data', data[200:210])
    introspect('uniqtext:', uniqtext)
    return [data,documents,sets,tags]


def do_this():
    check()
    infile = sys.argv[1]
    geneids = eval(open(sys.argv[2]).read())
    #print geneids
    exps_list = eval(open(sys.argv[3]).read())
    #size = sys.argv[2]
    #window = sys.argv[3]
    lines = open(infile).readlines()
    data = process_data(lines, geneids)
    documents = data[1]
    sets = data[2]
    tags = data[3]
    '''
    #print "data:", data
    print "num sets: ", len(sets.keys())
    print "sets:", sets
    print "sets[2]:", sets[2], len(sets[2])

    doc = sets[2][1]
    tagsfor2 = tags[2]
    print "doc:", doc
    print "tags for 2:", tagsfor2
    '''
    #print "documents: ", documents, "sets:", sets
    
    #test structure
    #print 'generif docs: ', documents[0:5] 
    
    
    #Train different models, by setting up different experimental settings here
    #exps_list = [{'size': 30, 'cw':2},{'size': 50, 'cw':2},{'size': 200, 'cw':2}]

    #Run deep learning training for each experimental setting
    for setting in exps_list:
        size = setting['size']
        window = setting['cw']
        orig_stdout = sys.stdout
        suffix = os.path.basename(infile)
        model_folder = 'model_gbg_vs%s_cw%s_%s' % (size,window,suffix.replace('.txt',''))
        os.system('mkdir %s' % model_folder)
        #create log file
        f = file('%s/model_gbg_vs%s_cw%s_%s.log.txt' % (model_folder,size,window,suffix), 'w')
        sys.stdout = f
        print "Model params: vector size - %s, window - %s, corpus - %s" % (size,window,suffix)
        print "length data: ", len(data), " length docs:", len(documents)
        print "length data[0]: ", len(data[0]), " length data[1]:", len(data[1])
        
        tic = time.time()
        for g in geneids:
            print "\nProcessing gene (ID): %s\n" % g
            docset = sets[g]
            tagset = tags[g]
            print 'Gene ID - "%s" (No. of RIFs - %s):' % (g, len(tagset))
            #print tagset
            model = doc2vec.Doc2Vec(docset, min_count = 1, window = int(window), size = int(size), workers=20, \
                                    dm_concat = 1,  dm_tag_count = 4 , dbow_words = 1) 
            #n = 1 #test
            #docvec = model.docvecs[n] 
            #print 'docvec  %s: ' % n , docvec, len(docvec)
            
            #Pick the last desc from set
            query_tag = tagset[-1][0] #[['SENT_11036822-9606-2'
            pmid = re.findall('([0-9]{8})', query_tag)[0]
            
            '''
            RESC02NW2PMG3QP:generank mxd074$ grep '9606\t635\t20628086' data/generif10-50K.txt 
            9606    635    20628086    2010-09-15 22:05:00    Observational study of gene-disease association, gene-environment interaction, and pharmacogenomic / toxicogenomic. (HuGE Navigator)
            RESC02NW2PMG3QP:generank mxd074$ grep '9606\t9\t17675654' data/generif10-50K.txt 
            9606    9    17675654    2008-03-13 09:01:00    Meta-analysis of gene-disease association and gene-environment interaction. (HuGE Navigator)
            9606    9    17675654    2010-01-21 00:00:00    systematic, literature-based review of the individual effects of NAT1 and NAT2 and their joint effects with smoking on bladder carcinogenesis
            '''
            os.system('grep "9606\t%s\t%s" data/generif10-50K.txt | cut -f 5 > tmp.txt' % (g,pmid))
            record = str(open('tmp.txt').read()).strip()
            print "\nQuery: Gene - %s, Pubmed ID - %s, RIF - %s" % (g, pmid, record)
           
            #Get similar records
            query_matches = model.docvecs.most_similar(query_tag)
            #print "Model similarities for query tag: %s" % query_tag, query_matches, len(model.docvecs.most_similar(query_tag))
            print "\nMost similar to query RIF - (%s, %s):" % (query_tag, record)
            for qm in query_matches:
                #print "Tag: ", qm[0]
                pmid = re.findall('([0-9]{8})', qm[0])[0]
                distance = qm[1]
                os.system('grep "9606\t%s\t%s" data/generif10-50K.txt | cut -f 5 > tmp.txt' % (g,pmid))
                record = str(open('tmp.txt').read()).strip()
                print '(PMID - %s, Distance - %s, RIF - %s)' % (pmid, distance, record)
            #create vector representation files
            modelfile = '%s/%s' % (model_folder,'model_gbg_vs%s_cw%s_g%s' % (size,window,g))
            write_model(model,modelfile,suffix,docset,tagset)
        toc = time.time()
        te = toc - tic 
        print "Time elapsed for training model: ", te 
        sys.stdout = orig_stdout
        f.close()

    return

do_this()
