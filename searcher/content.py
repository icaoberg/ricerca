# -*- coding: cp949 -*-

# Author: Ivan Cao-Berg, Baek Cho and Jennifer Bakal
# Created: December, 2011
#
# Copyright (C) 2012 Murphy Lab
# Lane Center for Computational Biology
# School of Computer Science
# Carnegie Mellon University
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

def distance( alpha, candidate, goodSet ):
 '''
 Calculates the distance between a candidate and every member
 of the good set.
 @return distance
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
 Calculate the norm between vector A and B.
 @param A
 @param B
 @alpha
 @return norm
 '''
 alpha = numpy.float64(1.0*alpha)
 A = numpy.float64( numpy.array( A ) )
 B = numpy.float64( numpy.array( B ) )
 return numpy.float64(numpy.power(numpy.sum(numpy.abs((A-B))**alpha),numpy.float64(1.0/alpha)))

def featnorm(trainset, testset):
 '''
 Feature normalization.
 @param train set
 @param test set
 @return normalized train and test sets
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
 @param train set
 @param test set
 @return normalized train and test sets
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
 Returns a ranked list.
 @param alpha
 @param candidates
 @param good set
 @return ranked list
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



