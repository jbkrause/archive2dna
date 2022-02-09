import array
from collections import defaultdict
from . import dna

class Representation:
    """Represents data array structured by nucleotides. Each column represents a DNA segment.
       Columns are indexed : 
           - they are dictionaries of two elements: the index and the actual data
           - columns are accessed using their index (not their position)
           - column_index containse the mapping between index and list postion
           - to move a column: change its index and reindex columns()"""
    def __init__(self,
                data_bytes = None,
                data_dna = None,
                n_lines   = 5,
                n_columns = 20 ):
        
        self.size = [n_lines, n_columns]
        self.column_index = {}
        self.column_index_min = 0
        self.column_index_max = n_columns
        
        # loading from bytes
        if data_bytes is not None:
        
            if len(data_bytes) > n_lines*n_columns:
                raise 

            # Data is organized by columns at initialization
            # each column is filled up sequentially
            # the last column is padded with zeroes at the end
            # (if not full)
            self.data = []
            for i in range(n_columns):
                i_from = i*n_lines
                i_to = (i+1)*n_lines    
                # if column will not be fully filled
                if i_to > len(data_bytes):
                    i_to = len(data_bytes)
                    delta = i_to - len(data_bytes)
                    for j in range(delta): # padding with zeros
                        self.data[-1]['column'].append(0)
                self.data.append( {'index':i,
                                   'column':array.array('b', 
                                       list( data_bytes[ i_from : i_to  ] ) )} )
            self.index_columns_num_currens()
                
        # loading from dna
        if data_dna is not None:
            # Data is organized by columns at initialization
            # each column corresponds to a DNA segment
            # - if a column of median size, it is padded using zeros
            # - if a column is longer thant median size it is imported as is
            self.data = []
            for i in range(n_columns):
                self.data.append( {'index':i, 'column':array.array('b') } )
                for j in range(len(data_dna[i])):
                    self.data[-1]['column'].append( dna.dna2bits( data_dna[i][j] ) )
                if len(data_dna[i]) < n_lines: # padding with zeroes
                    if i != n_columns-1: # not for last segement that is shorter
                        delta = n_lines - len(data_dna[i])
                        for j in range(delta):
                            self.data[-1]['column'].append(0)
            self.index_columns_num_currens()
            
    def index_columns_num_currens(self):
        """Indexes columns starting at 0 with increments of 1."""
        # This method is used for initial indexing
        self.column_index = defaultdict(int)
        for i in range( len(self.data) ):
            self.column_index[ self.data[i]['index'] ] = i
        self.column_index = dict(self.column_index)
        
    def reindex_columns(self):
        """Re-indexes columns, e.g. after loading DNA or inserting/removing a column."""
        self.column_index = {}
        for i in range( len(self.data) ):
            self.column_index[ self.data[i]['index'] ] = i
            
    def column_indexes(self):
        """Returns keys of column indexes, i.e. the actual column number that is
        used to acccess columns (not their internal position in representation)."""
        return self.column_index.keys()
                    
    def getcolumn(self, n, s=None):
        """Retruns whole column of index n."""
        if s==None:
            return self.data[ self.column_index[n]]['column']
        else:
            return self.data[ self.column_index[n]]['column'][s]
        
    def getpos(self, line, column):
        """Returns value at specific position in representation at specified line and column."""
        return self.data[ self.column_index[column]]['column'][line]

    #def setcolumn(self, n, column, start_at=0):
    #    for i in range(len(column)):
    #        self.data[ self.column_index[n]]['column'][i+start_at] = column[i]        

    def getline(self, n, s=None):
        line = array.array('b')
        if s == None:
            # TODO: performances, insure index sorted and do not sort here
            for i in sorted( self.column_indexes() ):
                col = self.getcolumn(i)
                if len(col) > n:
                    line.append( col[n] )
        else:
            raise "not implemented"
        return line

    #def setline(self, n, line, start_at=1):
    #    for i in range(len(line)):
    #        self.data[self.column_index[i+start_at]]['column'][n] = line[i]
            
    def setpos(self, line, column, value):
        """Sets value at specific position in representation at specified line and column."""
        self.data[ self.column_index[column]]['column'][line]=value

    def insertlines(self, position, n=1):
        """Inserts n lines at specified position."""
        for i in range(len(self.data)):
            col = self.data[i]['column']
            for j in range(n):
                col.insert(position,0)
        self.data[i]['column'] = col
        self.size[0] += n
        
    def insertcolumns(self, index, n=1):
        """Inserts n columns at specified index. 
        If any, exising indexes are shiftes"""
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
        
    def addcolumn(self, index):
        """Add a column at specified index. 
           - does NOT: check if columns already exist
           - does NOT: shift index of exisiting columns"""
        self.data.append({'index':index,
                          'column':array.array('b', [0]*self.size[0])})
        self.size[1] += 1
        self.reindex_columns()
        
    def popcolumn(self, index):
        """Removes column at index"""
        col = self.data.pop( self.column_index[index] )
        self.size[1] -= 1
        self.reindex_columns()
        return col
     
    def tonumpy(self):
        """
        Converts representation to numpy nd array.
        For debug purposes only. DO NOT USE IN LIBRARY."""
        import numpy as np
        out = np.array( np.full( self.size, None, dtype=object ) )
        for i in sorted( self.column_indexes() ):
            col = self.getcolumn(i)
            for j in range(len(col)):
                out[j,i] = col[j]
        return out
        
            
