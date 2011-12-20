import numpy
import scipy

def distance( alpha, candidate, goodSet, metric="norm" ):
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
  if metric == "norm":
   score = norm( candidate[2], goodSet[index][2] )

  if score != 0:
   score=numpy.float64(weight*numpy.exp(numpy.float64(alpha)*numpy.log(score)))
  elif candidate[1] < 0:
   flag_zero = True;
  else:
   score = 0
  
  total = total + score;
 
  #total=total[len(total)-1]
  #weights=weights[len(weights)-1]

 if flag_zero:
  total = numpy.float64(0)
 else:
  total = numpy.float64(total)/numpy.float64(weights)
  if total != 0:
   total=numpy.float64(numpy.exp(numpy.log(total)*(1.0/numpy.float64(candidate[1]))))
  elif numpy.float64(candidate[1])>0:
   total = 0
  else:
   if numpy.float64(candidate[1])<0:
    total = very_big

 return total 

def norm( A, B, alpha=2 ):
 alpha = numpy.float64(1.0*alpha)
 A = numpy.float64( numpy.array( A ) )
 B = numpy.float64( numpy.array( B ) )

 return numpy.float64(numpy.power(numpy.sum(numpy.abs((A-B))**alpha),numpy.float64(1.0/alpha)))

def ranking( alpha, candidates, goodSet ):

 #standard deviation of features in the dataset
 #std = numpy.std(numpy.array( goodSet[2]) )
 
 ratings = []
 iids = []
 for candidate in candidates:
  iids.append( candidate[0] )
  ratings.append( distance( alpha, candidate, goodSet ) )

 tups = zip(iids, ratings) # zip them as tuples

 result = sorted(tups, key=itemgetter(1))
 # note that this is sorting the tuples by the second element of the tuple
 # if you want to sort by the first element, then you should use ¡®itemgetter(0)¡¯
 sorted_iids = []
 sorted_scores = []
 for itm in result:
                sorted_iids.append(itm[0])
                sorted_scores.append(itm[1])

 return [sorted_iids, sorted_scores]


