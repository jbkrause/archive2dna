import array
from collections import defaultdict

class Representation:
    """Represents data array structured by nucleotides. Each column is a DNA segment.
       Columns are indexed : 
           - they are dictionaries of two elements: the index and the actual data
           - columns are accessed using their index (not their position in columns
             list via the column_index mapping between index and list postion
           - to move a column: change its index and update column_index"""
    def __init__(self,
                data_bytes = None,
                n_lines   = 5,
                n_columns = 20 ):
        
        self.size = [n_lines, n_columns]
        self.column_index = {}
        self.column_index_min = 0
        self.column_index_max = n_columns
        
        if len(data_bytes) > n_lines*n_columns:
            raise 

        # Data is organized by columns at initialization
        # each column is filled up sequentially
        self.data = []
        for i in range(n_columns):
            i_from = i*n_lines
            i_to = (i+1)*n_lines    
            # if column will not be fully filled
            if i_to > len(data_bytes):
                i_to = len(data_bytes)
            self.data.append( {'index':i,
                               'column':array.array('b', list( data_bytes[ i_from : i_to  ] ) )} )
            self.index_columns_num_currens()
            
    def index_columns_num_currens(self):
        self.column_index = defaultdict(int)
        for i in range( len(self.data) ):
            self.column_index[ self.data[i]['index'] ] = i
        self.column_index = dict(self.column_index)
        
    def reindex_columns(self):
        self.column_index = {}
        for i in range( len(self.data) ):
            self.column_index[ self.data[i]['index'] ] = i
                    
    def getcolumn(self, n, s=None):
        if s==None:
            return self.data[ self.column_index[n]]['column']
        else:
            return self.data[ self.column_index[n]]['column'][s]
        
    def getpos(self, line, column):
        return self.data[ self.column_index[column]]['column'][line]

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
            
    def setpos(self, line, column, value):
        self.data[ self.column_index[column]]['column'][line]=value
        
    def insertlines(self, position, n=1):
        """Inserts n lines at position"""
        for i in range(len(self.data)):
            col = self.data[i]['column']
            for j in range(n):
                col.insert(position,0)
        self.data[i]['column'] = col
        self.size[0] += n
        
    def insertcolumns(self, index, n=1):
        """Inserts n columns at index"""
        
        # shift indexes if required
        for i in range(len(self.data)):
            idx = self.data[i]['index']
            if idx >= index:
                self.data[i]['index'] += n
                
        # insert columns
        for i in range(index, index+n):
            self.data.append({'index':i,
                               'column':array.array('b', [0]*self.size[0])})
            
        self.size[1] += n
        self.reindex_columns()
