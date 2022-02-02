from unittest import TestCase
import os

import numpy as np

from archive2dna import dna
from archive2dna import package
from archive2dna import bytesutils

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
        binary_data = open('tests/data/aip_small.zip', 'rb').read()
        c = package.Container()
        c.load_binary(binary_data)
        self.assertTrue( isinstance(c.data, np.ndarray)   )

    def test_add_logical_redundancy(self):
        binary_data = open('tests/data/aip_small.zip', 'rb').read()
        c = package.Container()
        c.load_binary(binary_data)
        c.create_logical_redundancy()        
        self.assertTrue( isinstance(c.data, np.ndarray)   )
        
    def test_encode_write_decode_write(self):
        # encode bytes to dna and write
        binary_data = open('tests/data/aip_small.zip', 'rb').read()
        c = package.Container(package_id=None)
        c.load_binary(binary_data)
        c.add_outer_code()
        c.add_index()
        c.add_inner_code()
        c.to_dna()
        c.add_primers()
        text = c.write_dna()
        open('tests/tmp/dna_out.txt', 'w').write( text )
        c.compute_stats()
    
        # read dna, decode and write bytes
        c = package.Container()
        text = open('tests/tmp/dna_out.txt', 'r').read()
        c.read_dna(text)
        c.remove_primers()
        c.compute_segments_sizes()
        c.dna_to_array()
        c.dna_to_bits()
        c.decode_inner_code()
        c.sort_segments()
        c.decode_outer_code()
        binary_data = c.write_binary()
        open('tests/tmp/binary.out.zip', 'wb').write(binary_data)
        c.compute_stats()

        # check if input and output are the sha256 the same
        h1 = bytesutils.sha256('tests/data/aip_small.zip')
        h2 = bytesutils.sha256('tests/tmp/binary.out.zip')      
        # clean up and finish
        os.remove('tests/tmp/dna_out.txt')
        os.remove('tests/tmp/binary.out.zip')
        self.assertTrue( h1==h2  )


    def test_encode_write_decode_write_high_level(self): 
        # from bytes to DNA
        binary_data = open('tests/data/aip_small.zip', 'rb').read()
        c = package.Container(package_id='test:1')
        c.load_binary(binary_data) 
        c.create_logical_redundancy()
        c.convert_to_dna()
        text = c.write_dna()
        open('tests/tmp/dna_out.txt', 'w').write( text )
        c.compute_stats()

        # from DNA to bytes
        c = package.Container(package_id='test:1')
        text = open('tests/tmp/dna_out.txt', 'r').read()
        c.load_dna(text)
        c.check_and_correct_logical_redundancy()
        binary_data = c.write_binary()
        open('tests/tmp/binary.out.zip', 'wb').write(binary_data)
        c.compute_stats()

        # check if input and output are the sha256 the same
        h1 = bytesutils.sha256('tests/data/aip_small.zip')
        h2 = bytesutils.sha256('tests/tmp/binary.out.zip')     
        # clean up and finish
        os.remove('tests/tmp/dna_out.txt')
        os.remove('tests/tmp/binary.out.zip')
        self.assertTrue( h1==h2  )
        
