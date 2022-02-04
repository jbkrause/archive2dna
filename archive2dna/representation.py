import array
from collections import defaultdict

class Representation:

    def __init__(self,
                data_bytes = None,
                n_lines   = 5,
                n_columns = 20 ):
        
        self.size = [n_lines, n_columns]
        self.column_index = {}
        self.column_index_min = 1
        self.column_index_max = n_columns+1
        
        if len(data_bytes) > n_lines*n_columns:
            raise 

        self.data = []
        for i in range(n_columns):
            
            i_from = i*n_lines
            i_to = (i+1)*n_lines
            
            # if column will not be fully filled
            if i_to > len(data_bytes):
                i_to = len(data_bytes)
                
            self.data.append( {'index':i+1,
                               'column':array.array('b', list( data_bytes[ i_from : i_to  ] ) )} )
                             
            self.index_columns_num_currens()
            
    def index_columns_num_currens(self):
        self.column_index = defaultdict(int)
        for i in range( len(self.data) ):
            self.column_index[ self.data[i]['index'] ] = i
        self.column_index = dict(self.column_index)
                    
    def getcolumn(self, n, s=None):
        if s==None:
            return self.data[ self.column_index[n]]['column']
        else:
            return self.data[ self.column_index[n]]['column'][s]
        
    def setcolumn(self, n, column, start_at=0):
        for i in range(len(column)):
            self.data[ self.column_index[n]]['column'][i+start_at] = column[i]        

    def getline(self, n, s=None):
        line = array.array('b')
        if s == None:
            i_from = self.column_index_min
            i_to = self.column_index_max
        else:
            i_from = s.start
            i_to = s.stop
        for i in range( i_from, i_to ):
            if n < len( self.data[self.column_index[i]]['column'] ): # if column shorter
                line.append( self.data[self.column_index[i]]['column'][n] )
        return line

    def setline(self, n, line, start_at=1):
        for i in range(len(line)):
            self.data[self.column_index[i+start_at]]['column'][n] = line[i] 
            
    def insertcolumn(self, index, column):
        """Overwrites if index exists"""
        self.data.append({'index':index, 'column':column})
        self.column_index[index] = len(self.data)-1 
