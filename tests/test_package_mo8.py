from unittest import TestCase
import os

import numpy as np

from archive2dna import dna
from archive2dna import package
from archive2dna import bytesutils

# directories setup
test_package = 'tests/data/aip_olos.zip'
test_package = test_package.replace('/', os.sep)
test_tmp_dir = 'tests/tmp/'
test_tmp_dir = test_tmp_dir.replace('/', os.sep)
test_dna_tmp = test_tmp_dir + 'dna.txt'
test_aip_tmp = test_tmp_dir + 'aip.zip'
logging_file = test_tmp_dir + 'tests.log'

# saller outer code : test of blocks
mo = 8

if not os.path.isdir(test_tmp_dir):
    os.mkdir(test_tmp_dir)

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
        """Test mutation: replace bases in segments""" 
        # from bytes to DNA
        with open(test_package, 'rb') as f:
            binary_data = f.read()
        c = package.Container(package_id=None, mo=mo, logging_file=logging_file)
        c.load_binary(binary_data) 
        c.create_logical_redundancy()
        c.convert_to_dna()
        text = c.write_dna()
        with open(test_dna_tmp, 'w') as f:
            f.write( text )
            
        # do bases replacement
        with open(test_dna_tmp, 'r') as f:
            dna_segments = f.read().split('\n')
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
        
        with open(test_dna_tmp, 'w') as f:
            f.write( '\n'.join(dna_segments) )

        # from DNA to bytes
        c = package.Container(mo=mo, logging_file=logging_file)
        with open(test_dna_tmp, 'r') as f:
            text = f.read()
        c.load_dna(text)
        c.check_and_correct_logical_redundancy()
        binary_data = c.write_binary()
        with open(test_aip_tmp, 'wb') as f:
            f.write(binary_data)

        # check if input and output are the sha256 the same
        h1 = bytesutils.sha256(test_package)
        h2 = bytesutils.sha256(test_aip_tmp)             
        self.assertTrue( h1==h2  )
        
    def test_bases_deletion(self): 
        """Test mutation:  bases deletion in segments"""     
        # from bytes to DNA
        with open(test_package, 'rb') as f:
            binary_data = f.read()
        c = package.Container(package_id=None, mo=mo, logging_file=logging_file)
        c.load_binary(binary_data) 
        c.create_logical_redundancy()
        c.convert_to_dna()
        text = c.write_dna()
        with open(test_dna_tmp, 'w') as f:
            f.write( text )
        
        # do bases deletion
        with open(test_dna_tmp, 'r') as f:
            dna_segments = f.read().split('\n')
        # one
        dna_segments[200] = remove_base( dna_segments[200], pos=100 )
        with open(test_dna_tmp, 'w') as f:
            f.write( '\n'.join(dna_segments) )

        # from DNA to bytes
        c = package.Container(mo=mo, logging_file=logging_file)
        with open(test_dna_tmp, 'r') as f:
            text = f.read()
        c.load_dna(text)
        c.check_and_correct_logical_redundancy()
        binary_data = c.write_binary()
        with open(test_aip_tmp, 'wb') as f:
            f.write(binary_data)

        # check if input and output are the sha256 the same
        h1 = bytesutils.sha256(test_package)
        h2 = bytesutils.sha256(test_aip_tmp)
        self.assertTrue( h1==h2  )

    def test_bases_inseretion(self): 
        """Test mutation: insert bases in segments"""     
        # from bytes to DNA
        with open(test_package, 'rb') as f:
            binary_data = f.read()
        c = package.Container(package_id=None, mo=mo, logging_file=logging_file)
        c.load_binary(binary_data) 
        c.create_logical_redundancy()
        c.convert_to_dna()
        text = c.write_dna()
        with open(test_dna_tmp, 'w') as f:
            f.write( text )
        
        # do bases deletion
        with open(test_dna_tmp, 'r') as f:
            dna_segments = f.read().split('\n')
        # one
        dna_segments[200] = insert_base( dna_segments[250], pos=100 )
        with open(test_dna_tmp, 'w') as f:
            f.write( '\n'.join(dna_segments) )

        # from DNA to bytes
        c = package.Container(mo=mo, logging_file=logging_file)
        with open(test_dna_tmp, 'r') as f:
            text = f.read()
        c.load_dna(text)
        c.check_and_correct_logical_redundancy()
        binary_data = c.write_binary()
        with open(test_aip_tmp, 'wb') as f:
            f.write(binary_data)

        # check if input and output are the sha256 the same
        h1 = bytesutils.sha256(test_package)
        h2 = bytesutils.sha256(test_aip_tmp)
        self.assertTrue( h1==h2  )

    def test_segments_deletion(self): 
        """Test segments loss""" 
    
        # from bytes to DNA
        with open(test_package, 'rb') as f:
            binary_data = f.read()
        c = package.Container(package_id=None, mo=mo, logging_file=logging_file)
        c.load_binary(binary_data) 
        c.create_logical_redundancy()
        c.convert_to_dna()
        text = c.write_dna()
        with open(test_dna_tmp, 'w') as f:
            f.write( text )
        
        # do segments insertion
        with open(test_dna_tmp, 'r') as f:
            dna_segments = f.read().split('\n')
        for i in [10, 15, 24, 32, 51]:
            dna_sgements = remove_segments(dna_segments, i)
        with open(test_dna_tmp, 'w') as f:
            f.write( '\n'.join(dna_segments) )

        # from DNA to bytes
        c = package.Container(mo=mo, logging_file=logging_file)
        with open(test_dna_tmp, 'r') as f:
            text = f.read()
        c.load_dna(text)
        c.check_and_correct_logical_redundancy()
        binary_data = c.write_binary()
        with open(test_aip_tmp, 'wb') as f:
            f.write(binary_data)

        # check if input and output are the sha256 the same
        h1 = bytesutils.sha256(test_package)
        h2 = bytesutils.sha256(test_aip_tmp)             
        self.assertTrue( h1==h2  )


    def test_segments_permutations(self): 
        """Test segments permutations"""     
        # from bytes to DNA
        with open(test_package, 'rb') as f:
            binary_data = f.read()
        c = package.Container(package_id=None, mo=mo, logging_file=logging_file)
        c.load_binary(binary_data) 
        c.create_logical_redundancy()
        c.convert_to_dna()
        text = c.write_dna()
        with open(test_dna_tmp, 'w') as f:
            f.write( text )
        
        # do segments permutation
        with open(test_dna_tmp, 'r') as f:
            dna_segments = f.read().split('\n')
        dna_sgements = switch_segments(dna_segments, n1=10, n2=20)
        dna_sgements = switch_segments(dna_segments, n1=15, n2=8)
        dna_sgements = switch_segments(dna_segments, n1=31, n2=23)
        with open(test_dna_tmp, 'w') as f:
            f.write( '\n'.join(dna_segments) )

        # from DNA to bytes
        c = package.Container(mo=mo, logging_file=logging_file)
        with open(test_dna_tmp, 'r') as f:
            test = f.read()
        c.load_dna(text)
        c.check_and_correct_logical_redundancy()
        binary_data = c.write_binary()
        with open(test_aip_tmp, 'wb') as f:
            f.write(binary_data)

        # check if input and output are the sha256 the same
        h1 = bytesutils.sha256(test_package)
        h2 = bytesutils.sha256(test_aip_tmp)
        self.assertTrue( h1==h2  )

