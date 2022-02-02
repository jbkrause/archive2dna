#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This file is part of archive2dna.
#
# archive2dna is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Foobar is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Foobar. If not, see <https://www.gnu.org/licenses/>
#
# Author : Jan Krause-Bilvin
# First release: 2022-02-02

import os
import sys
import pprint

from archive2dna import package

usage = '''cli.py ACTION FILE_IN FILE_OUT [PACKAGE_ID]
           action : encode | decode\n'''

pp = pprint.PrettyPrinter(depth=6)

if len(sys.argv)==5:
    package_id = sys.argv[4]
    primer_length = 5
else:
    package_id = None
    primer_length = 0
    
if sys.argv[1] in ['-h', '--help']:
   print(usage)
    
elif sys.argv[1]=='encode':
    binary = sys.argv[2]
    dna = sys.argv[3]
    binary_data = open(binary, 'rb').read()
    c = package.Container( package_id=package_id, 
                           primer_length=primer_length,
                           necso=None )
    c.load_binary(binary_data) 
    c.create_logical_redundancy()
    c.convert_to_dna()
    c.compute_segments_sizes()
    text = c.write_dna()
    open(dna, 'w').write( text )
    pp.pprint(c.compute_stats())
    
elif sys.argv[1]=='decode':
    binary = sys.argv[3]
    dna = sys.argv[2]
    c = package.Container( package_id=package_id, 
                           primer_length=primer_length,
                           necso=None )
    #c.debug_output = True
    text = open(dna, 'r').read()
    c.load_dna(text)
    c.check_and_correct_logical_redundancy()
    binary_data = c.write_binary()
    open(binary, 'wb').write(binary_data)
    pp.pprint(c.compute_stats())

else:
    print(usage)


        
