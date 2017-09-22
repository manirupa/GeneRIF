#!/usr/local/bin/python

"""
@author: Manirupa Das
This script normalizes data, and given experimental settings, performs 
doc2vec training on sets of GenRIFs across dataset of genes (between-gene)

Additionally find the top-K, within-gene and across-gene similarities
for this type of training
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

# matlab integration
import matlab.engine
eng = matlab.engine.start_matlab()

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

def process_data(lines, origlines, geneids):
    data = []
    documents = []
    #also makes sets by geneid
    sets = {}
    tags = {}
    descs = {}
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
            origrecord = origlines[i]
            print record, origrecord
            label = 'SENT_%s-%s-%s' % (pubmedid,taxid, geneid)
            thewords = text.split()
            data.append([thewords,[label]])
            lbldsent = LabeledSentence(words=thewords, tags=[label,pubmedid,taxid,geneid])
            geneid = int(geneid)
            genefilepath = 'processed_genes/gene%s.txt' % geneid
            try:
               (sets[geneid]).append(lbldsent)            
            except:
                print "yoohoo"
                sets[geneid] = []
                (sets[geneid]).append(lbldsent)    
            try:
                (tags[geneid]).append([label,pubmedid,taxid,geneid])
            except:
                #print "yoohoo"
                tags[geneid] = []
                (tags[geneid]).append([label,pubmedid,taxid,geneid])
            try:
                temp = origrecord.strip().split('\t')
                desc = 'PMID - %s, Gene - %s, RIF - %s' % (temp[2], temp[1], temp[4])
                (descs[geneid]).append(desc)
            except:
                descs[geneid] = []
                temp = origrecord.strip().split('\t')
                desc = 'PMID - %s, Gene - %s, RIF - %s' % (temp[2], temp[1], temp[4])
                (descs[geneid]).append(desc)
            documents.append(lbldsent)
        else:
            pass
    #introspect('uniqtext:', uniqtext[:-10])
    return [data,documents,sets,tags,descs]


def do_this():
    check()
    infile = sys.argv[1]
    origfile=re.sub(r'\d+-gram.','',infile)
    geneids = eval(open(sys.argv[2]).read())
    #print geneids
    exps_list = eval(open(sys.argv[3]).read())
    #size = sys.argv[2]
    #window = sys.argv[3]
    lines = open(infile).readlines()
    origlines = open(origfile).readlines()
    data = process_data(lines, origlines, geneids)
    documents = data[1]
    sets = data[2] #sets of RIFS for each gene
    tags = data[3] #set of corresponding tags
    descs = data[4]
    
    introspect("data[0:10]", data[0][0:10])
    introspect("docs[0:10]", documents[0:10])
    print("sets[10]", sets[10])
    print("tags[10]", tags[10])
    
    print("descs[10]", descs[10])
    print(len(sets[10]), len(descs[10]), len(tags[10]))

    #Train different models, by setting up different experimental settings here
    #exps_list = [{'size': 30, 'cw':2},{'size': 50, 'cw':2},{'size': 200, 'cw':2}]
    #Run deep learning training for each experimental setting
    for setting in exps_list:
        size = int(setting['size'])
        window = int(setting['cw'])
        orig_stdout = sys.stdout
        suffix = os.path.basename(infile)
        model_folder = 'model_ovl_vs%s_cw%s_%s' % (size,window,suffix.replace('.txt',''))
        os.system('mkdir %s' % model_folder)
        #create log file
        f = file('%s/model_ovl_vs%s_cw%s_%s.log.txt' % (model_folder,size,window,suffix), 'w')
        sys.stdout = f
        print "Model params: vector size - %s, window - %s, corpus - %s" % (size,window,suffix)
        print "length data: ", len(data), " length docs:", len(documents)
        print "length data[0]: ", len(data[0]), " length data[1]:", len(data[1])
        
        tic = time.time()
        #train the model over all genes for this experiment setting
        model = doc2vec.Doc2Vec(documents, min_count = 1, window = window, size = size, workers=20, \
                                    dm_concat = 1,  dm_tag_count = 4 , dbow_words = 1) 
        toc = time.time()
        te = toc - tic 
        
        for g in geneids: #For each gene, produce model vector representation file for its docset, and similarities
            print "\nProcessing gene (ID): %s\n" % g
            docset = sets[g]
            tagset = tags[g]
            print 'Gene ID - "%s" (No. of RIFs - %s):' % (g, len(tagset))
            
            #Query GeneRIF-pick the last RIF desc from this gene set
            query_tag = tagset[-1][0] #[['SENT_11036822-9606-2'
            pmid = re.findall('([0-9]{8})', query_tag)[0]
            
            os.system('grep "9606\t%s\t%s" data/generif10-50K.txt | cut -f 5 > tmp.txt' % (g,pmid))
            record = str(open('tmp.txt').read()).strip()
            print "\nQuery: Gene - %s, Pubmed ID - %s, RIF - %s" % (g, pmid, record)
           
            #Get similar records - THIS RETRIEVES THE ACROSS-GENE SIMILARITIES, AS THIS IS ACROSS-GENE TRG 
            query_matches = model.docvecs.most_similar(query_tag)
            #print "Model similarities for query tag: %s" % query_tag, query_matches, len(model.docvecs.most_similar(query_tag))
            print "\nMost similar RIFs (BETWEEN-GENES) to query RIF - (%s, %s):\n" % (query_tag, record)
            for qm in query_matches:
                #print "Tag: ", qm[0]
                tag = qm[0]
                pmid = re.findall('([0-9]{8})', tag)[0] #extract the pubmed id from the tag
                distance = qm[1] #get the cosine distance from same tuple
                #try to get gene id
                try:
                    #tag = 'SENT_26217017-9606-2213'
                    #g = re.findall('(-[0-9]+$)', tag)[0]
                    #g.strip('-')
                    gene = re.findall('(-[0-9]+$)', tag)[0].strip('-')
                except: #use the first record
                    print "first record"
                    os.system('grep "%s" data/generif10-50K.txt | cut -f 2 | head -1 > tmp.txt' % pmid)
                    gene = str(open('tmp.txt').read()).strip()                    
                os.system('grep "9606\t%s\t%s" data/generif10-50K.txt | cut -f 5 | head -1 > tmp.txt' % (gene,pmid))
                record = str(open('tmp.txt').read()).strip()
                print '(PMID - %s, Distance - %s, Gene - %s, RIF - %s)' % (pmid, distance, gene, record)
            #create vector representation files
            modelfile = '%s/%s' % (model_folder,'model_ovl_vs%s_cw%s_g%s' % (size,window,g))
            print "\nModel file for settings (%s): %s" % (setting,modelfile)
            write_model(model,modelfile,suffix,docset,tagset)
            
            #NOW GET WITHIN-GENE SIMILARITIES FROM THIS MODEL FILE
            modelfile = '%s_%s' % (modelfile,suffix)            
            lines = open(modelfile).readlines() #Get vector reps for this gene
            query_idx = len(tagset) - 1
            print "Query GeneRIF index - %s" % query_idx
            query_tag = tags[g][query_idx][0]
            record = descs[g][query_idx]
            print "Query GeneRIF: %s \n" % record

            K = 10
            ret = eng.cosine(modelfile,query_idx,K) #Get top-K results
            indices = []
            distances = []
            #Now just get records from the indices
            for i in ret:
                rec = i.split("=")
                indices.append(int(rec[0]))
                distances.append(float(rec[1].strip()))
            print "Indices", indices[0:K]
            print "Distances", distances[0:K]
            
            topk_indices = indices[0:K]
            print "\nMost similar RIFs (WITHIN-GENE) to query RIF - (%s, %s):\n" % (query_tag, record)
            for i in range(K):
                print "(Distance: %s, Record: %s)" % (distances[i], descs[g][topk_indices[i]])

        print "Time elapsed for training model: ", te 
        sys.stdout = orig_stdout
        f.close()

    return

do_this()

eng.quit()
