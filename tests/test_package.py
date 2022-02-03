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

class PackageModudle(TestCase):


    def test_is_instanciated(self):
        c = package.Container() 
        self.assertTrue( isinstance(c.N, int) )

    def test_rand_mask(self):
        """Appling twice the same mask using XOR must result in identity"""
        c = package.Container() 
        b = b''
        for i in range(256):
            b += dna.int2bytes(i, n=1)        
        self.assertTrue( b == c.mask_bytes( c.mask_bytes(b) ) )
        
    def test_primers_management(self):
        """Tests identity after adding and removing primers"""
        sequence = 'A'
        c = package.Container(package_id='test:1')
        c.dna = [sequence]
        c.add_primers()
        c.remove_primers()
        self.assertTrue( c.dna[0] == sequence )

    def test_load_binary(self):
        with open(test_package, 'rb') as f:
            binary_data = f.read()
        c = package.Container()
        c.load_binary(binary_data)
        self.assertTrue( isinstance(c.data, np.ndarray)   )

    def test_add_logical_redundancy(self):
        with open(test_package, 'rb') as f:
            binary_data = f.read()
        c = package.Container()
        c.load_binary(binary_data)
        c.create_logical_redundancy()        
        self.assertTrue( isinstance(c.data, np.ndarray)   )
        
    def test_encode_write_decode_write(self):
        # encode bytes to dna and write
        with open(test_package, 'rb') as f:
            binary_data = f.read()
        c = package.Container(package_id=None)
        c.load_binary(binary_data)
        c.add_outer_code()
        c.add_index()
        c.add_inner_code()
        c.to_dna()
        c.add_primers()
        text = c.write_dna()
        with open(test_dna_tmp, 'w') as f:
            f.write( text )
        c.compute_stats()
    
        # read dna, decode and write bytes
        c = package.Container()
        with open(test_dna_tmp, 'r') as f:
            text = f.read()
        c.read_dna(text)
        c.remove_primers()
        c.compute_segments_sizes()
        c.dna_to_array()
        c.dna_to_bits()
        c.decode_inner_code()
        c.sort_segments()
        c.decode_outer_code()
        binary_data = c.write_binary()
        with open(test_aip_tmp, 'wb') as f:
            f.write(binary_data)
        c.compute_stats()

        # check if input and output are the sha256 the same
        h1 = bytesutils.sha256(test_package)
        h2 = bytesutils.sha256(test_aip_tmp)      
        self.assertTrue( h1==h2  )


    def test_encode_write_decode_write_high_level(self): 
        # from bytes to DNA
        with open(test_package, 'rb') as f:
            binary_data = f.read()
        c = package.Container(package_id='test:1')
        c.load_binary(binary_data) 
        c.create_logical_redundancy()
        c.convert_to_dna()
        text = c.write_dna()
        with open(test_dna_tmp, 'w') as f:
            f.write( text )
        c.compute_stats()

        # from DNA to bytes
        c = package.Container(package_id='test:1')
        with open(test_dna_tmp, 'r') as f:
            test = f.read()
        c.load_dna(text)
        c.check_and_correct_logical_redundancy()
        binary_data = c.write_binary()
        with open(test_aip_tmp, 'wb') as f:
            f.write(binary_data)
        c.compute_stats()

        # check if input and output are the sha256 the same
        h1 = bytesutils.sha256(test_package)
        h2 = bytesutils.sha256(test_aip_tmp)     
        self.assertTrue( h1==h2  )
        
