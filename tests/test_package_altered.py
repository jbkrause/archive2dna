from unittest import TestCase
import os

import numpy as np

from archive2dna import dna
from archive2dna import package
from archive2dna import bytesutils

def replace_base(s, pos=10, by='A'):
    return s[:pos] + by + s[(pos+1):]
    
def remove_base(s, pos=10):
    return s[:pos] + s[(pos+1):]

def insert_base(s, pos=10, by='A'):
    return s[:pos] + by + s[(pos):]

def switch_segments(segments, n1=10, n2=20):
    tmp = segments[n1]
    segments[n1]=segments[n2]
    segments[n2]=tmp
    return segments

def remove_segments(segments, n=100):
    segments = segments[:n] + segments[(n+1):]
    return segments

class PackageModudle(TestCase):

    def test_bases_replacement(self): 
    
        # from bytes to DNA
        binary_data = open('tests/data/aip_small.zip', 'rb').read()
        c = package.Container(package_id=None)
        c.load_binary(binary_data) 
        c.create_logical_redundancy()
        c.convert_to_dna()
        text = c.write_dna()
        open('tests/tmp/dna_out.txt', 'w').write( text )
        
        # do bases replacement
        dna_segments = open('tests/tmp/dna_out.txt', 'r').read().split('\n')
        # one replacement
        dna_segments[2] = replace_base( dna_segments[2], pos=2 )
        # one replacement
        dna_segments[10] = replace_base( dna_segments[10], pos=-4 )
        # two replacements congiguous
        dna_segments[20] = replace_base( dna_segments[20], pos=12 )
        dna_segments[20] = replace_base( dna_segments[20], pos=13 )
        # three replacements congiguous
        dna_segments[30] = replace_base( dna_segments[30], pos=20 )
        dna_segments[30] = replace_base( dna_segments[30], pos=21 )
        dna_segments[30] = replace_base( dna_segments[30], pos=22 ) 
        # four replacements congiguous
        dna_segments[35] = replace_base( dna_segments[35], pos=35 )
        dna_segments[35] = replace_base( dna_segments[35], pos=36 )
        dna_segments[35] = replace_base( dna_segments[35], pos=37 )
        dna_segments[35] = replace_base( dna_segments[35], pos=38 ) 
        # two replacements not-congiguous
        #dna_segments[40] = replace_base( dna_segments[40], pos=5 )
        #dna_segments[40] = replace_base( dna_segments[40], pos=18 )
        open('tests/tmp/dna_out.txt', 'w').write( '\n'.join(dna_segments) )

        # from DNA to bytes
        c = package.Container()
        text = open('tests/tmp/dna_out.txt', 'r').read()
        c.load_dna(text)
        c.check_and_correct_logical_redundancy()
        binary_data = c.write_binary()
        open('tests/tmp/binary.out.zip', 'wb').write(binary_data)

        # check if input and output are the sha256 the same
        h1 = bytesutils.sha256('tests/data/aip_small.zip')
        h2 = bytesutils.sha256('tests/tmp/binary.out.zip')
             
        # clean up and finish
        os.remove('tests/tmp/dna_out.txt')
        #os.remove('tests/tmp/binary.out.zip')
        self.assertTrue( h1==h2  )
        
    def test_bases_deletion(self): 
    
        # from bytes to DNA
        binary_data = open('tests/data/aip_small.zip', 'rb').read()
        c = package.Container(package_id=None)
        c.load_binary(binary_data) 
        c.create_logical_redundancy()
        c.convert_to_dna()
        text = c.write_dna()
        open('tests/tmp/dna_out.txt', 'w').write( text )
        
        # do bases deletion
        dna_segments = open('tests/tmp/dna_out.txt', 'r').read().split('\n')
        # one
        #dna_segments[2] = remove_base( dna_segments[2], pos=2 )
        # one
        #dna_segments[10] = remove_base( dna_segments[10], pos=-4 )
        # two  congiguous
        #dna_segments[20] = remove_base( dna_segments[20], pos=12 )
        #dna_segments[20] = remove_base( dna_segments[20], pos=12 )
        # three congiguous
        #dna_segments[30] = remove_base( dna_segments[30], pos=20 )
        #dna_segments[30] = remove_base( dna_segments[30], pos=20 )
        #dna_segments[30] = remove_base( dna_segments[30], pos=20 ) 
        # four congiguous
        #dna_segments[36] = remove_base( dna_segments[36], pos=35 )
        #dna_segments[36] = remove_base( dna_segments[36], pos=35 )
        #dna_segments[36] = remove_base( dna_segments[36], pos=35 )
        #dna_segments[36] = remove_base( dna_segments[36], pos=35 ) 
        # two not-congiguous
        #dna_segments[40] = remove_base( dna_segments[40], pos=5 )
        #dna_segments[40] = remove_base( dna_segments[40], pos=18 )
        open('tests/tmp/dna_out.txt', 'w').write( '\n'.join(dna_segments) )

        # from DNA to bytes
        c = package.Container()
        text = open('tests/tmp/dna_out.txt', 'r').read()
        c.load_dna(text)
        c.check_and_correct_logical_redundancy()
        binary_data = c.write_binary()
        open('tests/tmp/binary.out.zip', 'wb').write(binary_data)

        # check if input and output are the sha256 the same
        h1 = bytesutils.sha256('tests/data/aip_small.zip')
        h2 = bytesutils.sha256('tests/tmp/binary.out.zip')
             
        # clean up and finish
        os.remove('tests/tmp/dna_out.txt')
        #os.remove('tests/tmp/binary.out.zip')
        self.assertTrue( h1==h2  )


    def test_segments_deletion(self): 
    
        # from bytes to DNA
        binary_data = open('tests/data/aip_small.zip', 'rb').read()
        c = package.Container(package_id=None)
        c.load_binary(binary_data) 
        c.create_logical_redundancy()
        c.convert_to_dna()
        text = c.write_dna()
        open('tests/tmp/dna_out.txt', 'w').write( text )
        
        # do segments deletion
        dna_segments = open('tests/tmp/dna_out.txt', 'r').read().split('\n')
        for i in [10, 15, 24, 32, 51]:
            dna_sgements = remove_segments(dna_segments, i)
        open('tests/tmp/dna_out.txt', 'w').write( '\n'.join(dna_segments) )

        # from DNA to bytes
        c = package.Container()
        text = open('tests/tmp/dna_out.txt', 'r').read()
        c.load_dna(text)
        c.check_and_correct_logical_redundancy()
        binary_data = c.write_binary()
        open('tests/tmp/binary.out.zip', 'wb').write(binary_data)

        # check if input and output are the sha256 the same
        h1 = bytesutils.sha256('tests/data/aip_small.zip')
        h2 = bytesutils.sha256('tests/tmp/binary.out.zip')
             
        # clean up and finish
        os.remove('tests/tmp/dna_out.txt')
        #os.remove('tests/tmp/binary.out.zip')
        self.assertTrue( h1==h2  )


    def test_segments_permutations(self): 
    
        # from bytes to DNA
        binary_data = open('tests/data/aip_small.zip', 'rb').read()
        c = package.Container(package_id=None)
        c.load_binary(binary_data) 
        c.create_logical_redundancy()
        c.convert_to_dna()
        text = c.write_dna()
        open('tests/tmp/dna_out.txt', 'w').write( text )
        
        # do segments permutations
        dna_segments = open('tests/tmp/dna_out.txt', 'r').read().split('\n')
        dna_sgements = switch_segments(dna_segments, n1=10, n2=20)
        dna_sgements = switch_segments(dna_segments, n1=15, n2=8)
        dna_sgements = switch_segments(dna_segments, n1=31, n2=23)
        open('tests/tmp/dna_out.txt', 'w').write( '\n'.join(dna_segments) )

        # from DNA to bytes
        c = package.Container()
        text = open('tests/tmp/dna_out.txt', 'r').read()
        c.load_dna(text)
        c.check_and_correct_logical_redundancy()
        binary_data = c.write_binary()
        open('tests/tmp/binary.out.zip', 'wb').write(binary_data)

        # check if input and output are the sha256 the same
        h1 = bytesutils.sha256('tests/data/aip_small.zip')
        h2 = bytesutils.sha256('tests/tmp/binary.out.zip')
             
        # clean up and finish
        os.remove('tests/tmp/dna_out.txt')
        #os.remove('tests/tmp/binary.out.zip')
        self.assertTrue( h1==h2  )

