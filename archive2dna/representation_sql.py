import array
from collections import defaultdict
from . import dna

from sqlalchemy import create_engine
from sqlalchemy import MetaData, Table, Column
from sqlalchemy import Integer, SmallInteger, String, and_, or_
from sqlalchemy import select, insert, update, delete


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
                numblocks = 1,
                dblocksize = 10,
                dnecso = 3,
                dN = 10,
                n_lines   = 5,
                n_columns = 20 ):
        
        self.size = [n_lines, n_columns]
        self.column_index = {}
        self.column_index_min = 0
        self.column_index_max = n_columns
        
        self.meta = MetaData()
        
        # Prepare spaceace for inner code and index
        delta_lines   = dN - n_lines
        self.size[0] += delta_lines

        args = ['representation', self.meta, Column('id', Integer, primary_key = True),  Column('index', Integer)]
        for i in range(delta_lines + n_lines):
            args.append( Column('c'+str(i), SmallInteger) )    
        representation = Table( *args )

        self.engine = create_engine('sqlite://', echo=False)
        self.meta.create_all(self.engine)
        self.table = self.meta.tables['representation']
     
        # loading from bytes
        if data_bytes is not None:
        
            if len(data_bytes) > n_lines*n_columns:
                raise 
                
            # Data is organized by columns at initialization
            # each column is filled up sequentially
            # the last column is padded with zeroes at the end
            # (if not full)
            # * dnecso comlumns are reserved for outer code
            # * delta_lines = dN-Nlines lines are reserved for inner code and index           

            idx = -1
            for blk in range(numblocks):
                for i in range(dnecso):
                    idx += 1
                    args = {'index':idx}
                    for ix, x in enumerate( [0 for i in range(dN)] ):
                        args[ 'c'+str(ix) ] = x
                    stmt = insert(self.table).values(**args)
                    with self.engine.connect() as connection:
                        connection.execute(stmt)                    
                for i in range(dblocksize - dnecso):
                    idx += 1
                    i_from = i*n_lines + (dblocksize-dnecso)*blk
                    i_to = (i+1)*n_lines + (dblocksize-dnecso)*blk   
                    # if column will not be fully filled
                    if i_to > len(data_bytes):
                        i_to = len(data_bytes)
                        delta = i_to - len(data_bytes)
                        for j in range(delta): # padding with zeros
                            self.data[-1]['column'].append(0)
                    args = {'index':idx}
                    for ix, x in enumerate( [0 for i in range(delta_lines)] + list( data_bytes[ i_from : i_to  ] )  ):
                        args[ 'c'+str(ix) ] = x
                    stmt = insert(self.table).values(**args)
                    with self.engine.connect() as connection:
                        connection.execute(stmt)
                        

#                for i in range(n_columns):
#                    idx += 1
#                    i_from = i*n_lines
#                    i_to = (i+1)*n_lines    
#                    # if column will not be fully filled
#                    if i_to > len(data_bytes):
#                        i_to = len(data_bytes)
#                        delta = i_to - len(data_bytes)
#                        for j in range(delta): # padding with zeros
#                            self.data[-1]['column'].append(0)
#                    args = {'index':i}
#                    for ix, x in enumerate( [None for i in range(delta_lines)] + list( data_bytes[ i_from : i_to  ] ) ):
#                        args[ 'c'+str(ix) ] = x
#                    stmt = insert(self.table).values(**args)
#                    with self.engine.connect() as connection:
#                        connection.execute(stmt)                
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
        
    #def reindex_columns(self):
    #    """Re-indexes columns, e.g. after loading DNA or inserting/removing a column."""
    #    self.column_index = {}
    #    for i in range( len(self.data) ):
    #        self.column_index[ self.data[i]['index'] ] = i
            
    def column_indexes(self):
        """Returns keys of column indexes, i.e. the actual column number that is
        used to acccess columns (not their internal position in representation)."""
        stmt = select([
            self.table.columns[ 'index' ] ]
            ).order_by( self.table.columns.index )
        with self.engine.connect() as connection:
            results = connection.execute(stmt).fetchall()
        return [x[0] for x in results]
                    
    def getcolumn(self, n, s=None):
        """Retruns whole column of index n.
        An optional slice s may be specified to restrict returned range."""
        #fields = ['c'+str(i) for i in range(self.dN)]
        stmt = select(
            self.table.columns
            ).where(and_(
                self.table.columns.index == n   
            ))
        with self.engine.connect() as connection:
            results = connection.execute(stmt).fetchall()
        if s == None:
            return results[0][2:]
        else:
            return results[0][2:][s]
        # FIXME: return array 'b' instead?
        
    def getpos(self, line, column):
        """Returns value at specific position in representation at specified line and column."""
        line_label = 'c'+str(line)
        stmt = select([
            self.table.columns[ line_label ] ]
            ).where(and_(
                self.table.columns.index == column
            ))             
        with self.engine.connect() as connection:
            res = connection.execute(stmt).fetchall()
        return [x[0] for x in res][0]
        #return self.data[ self.column_index[column]]['column'][line]   

    def getline(self, n, s=None):
        """Returns whole line n by default. 
        An optional slice s may be specified to restrict returned range."""
        if s == None:
            stmt = select([
                self.table.columns[ 'c'+str(n) ] ]
                ).order_by( self.table.columns.index )
        else:
            stmt = select([
                self.table.columns[ 'c'+str(n) ] ]
                ).where(and_(
                    self.table.columns.index >= s.start,
                    self.table.columns.index < s.stop
                )).order_by( self.table.columns.index )            
        with self.engine.connect() as connection:
            line = connection.execute(stmt).fetchall()
            #print(line)
        return [x[0] for x in line]
            
    def setpos(self, line, column, value):
        """Sets value at specific position in representation at specified line and column."""
        line_label = 'c'+str(line)
        stmt = update(self.table, values = { line_label : value } 
                ).where(and_(
                    self.table.columns.index == column
                ))         
        with self.engine.connect() as connection:
            connection.execute(stmt)
        #self.data[ self.column_index[column]]['column'][line]=value
        
    #def insertcolumns(self, index, n=1):
    #    """Inserts n columns at specified index. 
    #    If any, exising indexes are shiftes"""
    #    # shift indexes if required
    #    for i in range(len(self.data)):
    #        idx = self.data[i]['index']
    #        if idx >= index:
    #            self.data[i]['index'] += n         
    #    # insert columns
    #    for i in range(index, index+n):
    #        self.data.append({'index':i,
    #                           'column':array.array('b', [0]*self.size[0])})    
    #    self.size[1] += n
    #    self.reindex_columns()
        
    def addcolumn(self, index):
        """Add a column at specified index. 
           - does NOT: check if columns already exist
           - does NOT: shift index of exisiting columns"""
        args = {'index':index}
        for ix, x in enumerate( [0 for i in range(dN)] ):
            args[ 'c'+str(ix) ] = x
        stmt = insert(self.table).values(**args)
        with self.engine.connect() as connection:
            connection.execute(stmt)     
        #self.data.append({'index':index,
        #                  'column':array.array('b', [0]*self.size[0])})
        self.size[1] += 1
        #self.reindex_columns()
        
    def popcolumn(self, index):
        """Removes column at index"""
        stmt = delete(self.table 
                ).where(and_(
                    self.table.columns.index == index
                ))         
        with self.engine.connect() as connection:
            connection.execute(stmt)
        #col = self.data.pop( self.column_index[index] )
        self.size[1] -= 1
        #self.reindex_columns()
        return None
     
    def tonumpy(self):
        """
        Converts representation to numpy nd array.
        For debug purposes only. DO NOT USE IN LIBRARY."""
        import numpy as np
        stmt = select(self.table.columns).order_by(self.table.columns.index)
        with self.engine.connect() as connection:
            cols = connection.execute(stmt).fetchall()   
        out = np.array( np.full( [self.size[0], len(cols)], None, dtype=object ) ) # FIXME sietz[0] ?
        for icol, col in enumerate(cols):
            col2 = col[2:]#remove id and index
            for j in range(len(col2)): #range(len(col2)) range(self.size[0]) # FIXME sietz[0] ?
                out[j,icol] = col2[j]            
        return out
        
            
