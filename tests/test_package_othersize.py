from unittest import TestCase
import os

import numpy as np

from archive2dna import dna
from archive2dna import package
from archive2dna import bytesutils

# directories setup

test_tmp_dir = 'tests/tmp/'
test_tmp_dir = test_tmp_dir.replace('/', os.sep)
test_dna_tmp = test_tmp_dir + 'dna.txt'
test_aip_tmp = test_tmp_dir + 'aip.zip'

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

    def test_olos2(self): 
        """Test mutation: replace bases in segments""" 
        
        test_package = 'tests/data/aip_olos2.zip'
        test_package = test_package.replace('/', os.sep)
        
        # from bytes to DNA
        with open(test_package, 'rb') as f:
            binary_data = f.read()
        c = package.Container(package_id=None)
        c.load_binary(binary_data) 
        c.create_logical_redundancy()
        c.convert_to_dna()
        text = c.write_dna()
        with open(test_dna_tmp, 'w') as f:
            f.write( text )

        # from DNA to bytes
        c = package.Container()
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


    def test_matterhorn(self): 
        """Test mutation: replace bases in segments""" 
        
        test_package = 'tests/data/aip_matterhorn.zip'
        test_package = test_package.replace('/', os.sep)
        
        # from bytes to DNA
        with open(test_package, 'rb') as f:
            binary_data = f.read()
        c = package.Container(package_id=None)
        c.load_binary(binary_data) 
        c.create_logical_redundancy()
        c.convert_to_dna()
        text = c.write_dna()
        with open(test_dna_tmp, 'w') as f:
            f.write( text )

        # from DNA to bytes
        c = package.Container()
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

