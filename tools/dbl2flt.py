#! /usr/bin/env python3
'''
Converts VLSV file in-place from 8 to 4 byte floats.

Copyright 2026 Finnish Meteorological Institute

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


Author(s): Ilja Honkonen
'''

from sys import argv

from numpy import fromfile, zeros

with open(argv[1], 'rb') as infile:
	infile.seek(8)
	footer_start = fromfile(infile, dtype = 'uint64', count = 1)[0]
	infile.seek(footer_start)
	footer = infile.read().decode('ascii')
	lines = footer.split('\n')

with open(argv[1], 'r+b') as outfile:
	for line in lines:
		if not 'datatype="float"' in line:
			continue
		if not 'datasize="8"' in line:
			continue
		data_start = int(line.split('>')[1].split('<')[0])
		arr_size = int(line.split('arraysize=')[1].split('"')[1])
		vec_size = int(line.split('vectorsize=')[1].split('"')[1])
		data_size = arr_size * vec_size * 8
		for i in range(arr_size): # don't use much memory
			outfile.seek(data_start + 8*i*vec_size)
			data = fromfile(outfile, dtype = 'float64', count = vec_size)
			outfile.seek(data_start + 4*i*vec_size)
			data.astype('float32').tofile(outfile)
		outfile.seek(data_start + data_size//2)
		zeros(data_size//8).astype('uint32').tofile(outfile)
	outfile.seek(footer_start)
	for line in lines:
		if 'datatype="float"' in line and 'datasize="8"' in line:
			line = line.replace('datasize="8"', 'datasize="4"')
		outfile.write((line + '\n').encode('utf-8'))
