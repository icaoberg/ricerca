import cPickle as pickle

file = open( 'slf33.pkl' )
slf33 = pickle.load( file )
file.close()

database = []
for i in range(len(slf33[2])):
 database.append( [slf33[0][i], slf33[1][i], slf33[2][i]] )
 
file = open( 'database.pkl', 'w' )
slf33 = pickle.dump( database, file )
file.close()
