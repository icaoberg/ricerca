# -*- coding: cp949 -*-

# Author: Ivan Cao-Berg, Baek Cho and Jennifer Bakal
# Created: December, 2011
#
# Copyright (C) 2011-2013 Murphy Lab
# Lane Center for Computational Biology
# School of Computer Science
# Carnegie Mellon University
#
# February 27,2013 J. Bakal Moved combinePosandNeg function from
#                           searchContent.py/omero.searcher.py
#                           Added rankingWrapper
#
# March 27, 2013   J. Bakal Updated rankingWrapper
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation; either version 2 of the License,
# or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA.
#
# For additional information visit http://murphylab.web.cmu.edu or
# send email to murphy@cmu.edu

import numpy
import scipy
from operator import itemgetter, attrgetter

def loadContentDB(filename):
    '''
    Helper method that loads a content database from file on disk

    :param filename: the filename of a file containing valid content database
    :type alpha: string
    :rtype: content database
    '''

    try:
        import cPickle as pickle
    except:
        import pickle

    try:
        cdb = pickle.load(open(filename))
        return cdb
    except:
        print "Unable to load file: " + filename
        return {}

def getDBscales(cdb,query_scale):
    '''
    Helper method that returns the scales from a content database

    :param cdb: content database
    :type cdb: dictionary
    :rtype: scales
    '''

    keys=cdb.keys()
    keys.remove('info')
    dbkeys=[v for v in keys if v>.75*query_scale and v<1.5*query_scale]
    return dbkeys

def rankingWrapper(contentDB, image_refs_dict, processIDs, processSearchSet):
    '''
    Wrapper method that performs a ranking given a content database
    '''

    '''
    @param contentDB
    @param image_refs_dict
    @param processIDs
    @param processSearchSet
    @return final_result, dscale

    image_refs_dict has format image_refs_dict[image_ID]=[(image_scale,image_dna_reference_ID),similarity] 
        where image_ID is platform specific and similarity is 1 for images most similar and -1 for images most dissimilar

    processIDs is a function which takes a row from the contentDB and converts it to a platform specific ID to be used for output

    processSearchSet is a function which takes as input the contentDB, an image_refs_dict with format 
        image_refs_dict[image_id]=(image_scale,image_dna_reference_ID), and the contentDB scale determined from
        the scales of the images and returns a list formatted for the ranking search with each row of the format:
        [image_info_for_output,1,[features]]
    '''

    def divideImageDict(image_refs_dict):
        pos_image_dict={}
        neg_image_dict={}

        for key,val in image_refs_dict.items():
            if val[1]==1:
                pos_image_dict[key]=val[0]
            elif val[1]==-1:
                 neg_image_dict[key]=val[0]

        return pos_image_dict,neg_image_dict

    def rankSearchSet(contentDB,image_refs_dict,process_IDs,processSearchSet):
        #determine resolution to use                                                        
        keys = contentDB.keys()
        keys.remove('info')
        keys.sort()
        dscale = keys[0]

        scale = max([val[0] for val in image_refs_dict.values()])
        print 'scale from query images is',scale

        #find scale in the dictionary closest to the scale of the local images              
        for key in keys:
            if abs(key - scale) < abs(dscale - scale):
                dscale = key

        print 'scale/s from database is', dscale

        if not contentDB.has_key( dscale ):
            sys.exit("System error - scale not found in content database.")
        else:
            data=contentDB[dscale]

        dataset=[]
        for cdb_row in data:
            ID=processIDs(cdb_row)
            feat_vec = list(cdb_row[11:])
            dataset.append([ID,0,feat_vec])

        search_set=processSearchSet(contentDB, image_refs_dict, dscale)

        alpha = -5

        #RANKING IMAGES USING RICERCA                                                       
        print "Ranking images"
        normalization = 'zscore'
        if search_set:
            [sorted_iids,sorted_scores] = ranking( alpha, dataset, search_set, normalization )
            return [sorted_iids, sorted_scores, dscale]
        else:
             return [[],[],0]


    def combinePosandNeg(pos_sorted,neg_sorted):
        #reverse negative list before averaging                                             
        neg_sorted.reverse()

        img_avg=[]
        rank=range(len(pos_sorted))
        pos_sort_rank=zip(pos_sorted,rank)
        neg_sort_rank=zip(neg_sorted,rank)
        pos_sort_rank.sort(key=lambda img:img[0])
        neg_sort_rank.sort(key=lambda img:img[0])
        for j in range(len(pos_sort_rank)):
            if pos_sort_rank[j][0]==neg_sort_rank[j][0]:
                img_rank=(pos_sort_rank[j][0],float(pos_sort_rank[j][1]+neg_sort_rank[j][1])/2)
                img_avg.append(img_rank)

        #sort by rank                                                                       
        img_avg.sort(key=lambda img:img[1])

        avg_sorted=[]
        for img in img_avg:
            avg_sorted.append(img[0])

        return avg_sorted

    #separate positive and negative images                                                  
    pos_image_dict,neg_image_dict=divideImageDict(image_refs_dict)

    #rank database images based on positive images                                          
    if pos_image_dict:
        try:                                                                               
            [sorted_iids_pos, sorted_scores_pos, dscale] = rankSearchSet(contentDB, pos_image_dict, processIDs, processSearchSet)
        except:                                                                            
            return [],0                                                                   

    #rank database images based on negative images
    if neg_image_dict:
        try:
            [sorted_iids_neg, sorted_scores_neg, dscale] = rankSearchSet(contentDB, neg_image_dict, processIDs, processSearchSet)
        except:
            return [],0

    #combine &/or report rankings
    if pos_image_dict:
        if neg_image_dict: #there are both positive and negative samples
            final_result = combinePosandNeg(sorted_iids_pos,sorted_iids_neg)
            final_result = (final_result, None)
        else:              #there are only positive samples
            final_result = (sorted_iids_pos, sorted_scores_pos)
    elif neg_image_dict:     #there are only negative samples
        final_result = (sorted_iids_neg, sorted_scores_neg)

    return final_result,dscale


def distance( alpha, candidate, goodSet ):
 '''
 Calculates the distance between a candidate and every member of the good set

 :param alpha: alpha
 :type alpha: double
 :param candidate: a feature vector 
 :type candidates: list
 :param goodSet: a list of feature vectors
 :type goodSet: array
 :rtype: distance
 '''

 very_big = float(numpy.finfo( numpy.float32 ).max)/2;
 flag_zero = False
 total = numpy.float64(0)

 #number of images
 counts = len( goodSet )

 #pairwise distance
 weights = numpy.float64(0)

 for index in range( counts ):
  weight = numpy.float64(goodSet[index][1])
  weights = weights + weight
  d = norm( candidate[2], goodSet[index][2] )
  score = weight * numpy.power(d, numpy.float64(alpha)) 
  total = total + score;

 if candidate[1] < 0:
  total = numpy.float64(0)
 else:
  total = numpy.float64(total)/numpy.float64(weights)
  if total != 0:
      total=numpy.float64(numpy.power(total,numpy.float64(alpha)))
  elif numpy.float64(candidate[1])>0:
   total = 0
  else:
   if numpy.float64(candidate[1])<0:
    total = very_big

 return total 

def norm( A, B, alpha=2 ):
 '''
 Calculates the norm between vectors A and B

 :param A: a vector
 :type A: list of doubles
 :param B: a vector 
 :type B: list of doubles
 :rtype: norm between vectors A and B
 '''

 alpha = numpy.float64(1.0*alpha)
 A = numpy.float64( numpy.array( A ) )
 B = numpy.float64( numpy.array( B ) )
 return numpy.float64(numpy.power(numpy.sum(numpy.abs((A-B))**alpha),numpy.float64(1.0/alpha)))

def featnorm(trainset, testset):
 '''
 Feature normalization.

 :param trainset: training set
 :type trainset: list of feature vectors
 :param testset: test set 
 :type testset: list of feature vectors
 :rtype: normalized train and test sets
 '''

 trainset_id = []
 trainset_wt = []
 trainset_feat = []

 for itm in trainset:
  trainset_id.append(itm[0])
  trainset_wt.append(itm[1])
  trainset_feat.append(itm[2])

 trainset_feat = numpy.array(trainset_feat)
 min_col = trainset_feat.min(axis=0) + 1e-10
 max_col = trainset_feat.max(axis=0) + 1e-10

 trainset_normfeat = (trainset_feat-min_col)/numpy.float64(max_col) 
    
 testset_id = []
 testset_wt = []
 testset_feat = []

 for itm in testset:
  testset_id.append(itm[0])
  testset_wt.append(itm[1])
  testset_feat.append(itm[2])

 testset_feat = numpy.array(testset_feat)
 testset_normfeat = (testset_feat-min_col)/numpy.float64(max_col)

 new_trainset = []
 for i in range(len(trainset)):
  new_trainset.append([trainset_id[i], trainset_wt[i], trainset_normfeat[i]])

 new_testset = []
 for i in range(len(testset)):
  new_testset.append([testset_id[i], testset_wt[i], testset_normfeat[i]])

 return new_trainset, new_testset

def featnorm_z(trainset, testset):
 '''
 z-Score feature normalization.

 :param trainset: training set
 :type trainset: list of feature vectors
 :param testset: test set 
 :type testset: list of feature vectors
 :rtype: normalized train and test sets
 '''

 trainset_id = []
 trainset_wt = []
 trainset_feat = []

 for itm in trainset:
  trainset_id.append(itm[0])
  trainset_wt.append(itm[1])
  trainset_feat.append(itm[2])

 trainset_feat = numpy.array(trainset_feat)

 mean_col = trainset_feat.mean(axis=0)
 std_col = trainset_feat.std(axis=0) + 1e-10

 trainset_normfeat = (trainset_feat - mean_col)/numpy.float64(std_col)

 testset_id = []
 testset_wt = []
 testset_feat = []

 for itm in testset:
     testset_id.append(itm[0])
     testset_wt.append(itm[1])
     testset_feat.append(itm[2])

 testset_feat = numpy.array(testset_feat)
 testset_normfeat = (testset_feat-mean_col)/numpy.float64(std_col)

 new_trainset = []
 for i in range(len(trainset)):
     new_trainset.append([trainset_id[i], trainset_wt[i], trainset_normfeat[i]])

 new_testset = []
 for i in range(len(testset)):
     new_testset.append([testset_id[i], testset_wt[i], testset_normfeat[i]])

 return new_trainset, new_testset

def ranking( alpha, candidates, goodSet, normalization='zscore' ):
 '''
 Returns a ranked list

 :param alpha: alpha
 :type alpha: double
 :param candidates: list of image ids from candidates
 :type candidates: double
 :param goodSet: list of image ids from members of the good set
 :type goodSet: list of longs
 :param normalization: normalization parameter. default value is 'zscore'
 :type normalization: string
 :rtype: list of ranked image 
 '''

 #standard deviation of features in the dataset
 #std = numpy.std(numpy.array( goodSet[2]) )

 #normalize the feature vector
 if normalization == 'zscore':
  [candidates, goodSet] = featnorm_z(candidates, goodSet)
 else:
  [candidates, goodSet] = featnorm(candidates, goodSet)

 ratings = []
 iids = []
 for candidate in candidates:
  iids.append( candidate[0] )
  ratings.append( distance( alpha, candidate, goodSet ) )

 tups = zip(iids, ratings) # zip them as tuples
 
 result = sorted(tups, key=itemgetter(1))
 # note that this is sorting the tuples by the second element of the tuple
 # if you want to sort by the first element, then you should use ‘itemgetter(0)’
 sorted_iids = []
 sorted_scores = []
 for itm in result:
  sorted_iids.append(itm[0])
  sorted_scores.append(itm[1])

 return [sorted_iids, sorted_scores]



