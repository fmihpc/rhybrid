#! /usr/bin/env python3
'''
Sorts RHybrid VLSV cell data in-place into ascending cell id order.

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

from xml.parsers import expat
from numpy import fromfile

current_name = None

# Returns dict of variable,etc info and data offsets
def load_metadata(filename):
	sim_data = dict()
	sim_data['filename'] = filename

	# read footer
	with open(filename, 'rb') as infile:
		infile.seek(8)
		footer_start = fromfile(infile, dtype = 'uint64', count = 1)[0]
		infile.seek(footer_start)
		footer = infile.read().decode('ascii')


	def element_start(name, attrs):
		global current_name
		if name == 'VARIABLE' \
			or name == 'PARAMETER' \
			or name.startswith('MESH_NODE_CRDS_') \
		:
			if current_name != None:
				exit('internal error in element_start()')
			if not name.startswith('MESH_NODE_CRDS_'):
				name = attrs['name']
			current_name = name
			sim_data[name] = attrs

	def element_end(name):
		global current_name
		if name == 'VARIABLE' \
			or name == 'PARAMETER' \
			or name.startswith('MESH_NODE_CRDS_') \
		:
			if current_name == None:
				exit('internal error in element_end()')
			current_name = None

	def char_data(data):
		global current_name
		if data.strip() != '':
			if current_name != None:
				sim_data[current_name]['fileoffset'] = int(data)

	parser = expat.ParserCreate()
	parser.StartElementHandler = element_start
	parser.EndElementHandler = element_end
	parser.CharacterDataHandler = char_data
	parser.Parse(footer, True)
	return sim_data


# Adds to given dict array data of variables,etc in given dict
def load_data(data):
	with open(data['filename'], 'rb') as infile:
		for key in data.keys():
			if not 'type' in data[key]:
				continue
			if not data[key]['type'] == 'celldata':
				continue

			dtype = None
			if data[key]['datatype'] == 'float':
				if data[key]['datasize'] == '8':
					dtype = 'double'
				else:
					exit('float size != 8 not supported')
			elif data[key]['datatype'] == 'uint':
				if data[key]['datasize'] == '8':
					dtype = 'uint64'
				elif data[key]['datasize'] == '4':
					dtype = 'uintc'
				else:
					exit('float size != 8 not supported')
			if dtype == None:
				exit('unsupported datatype: ' + data[key]['datatype'])
			if data[key]['vectorsize'] != '1':
				dtype = data[key]['vectorsize'] + dtype

			infile.seek(data[key]['fileoffset'])
			data[key]['data'] = fromfile(
				infile, dtype = dtype, count = int(data[key]['arraysize']))
			if len(data[key]['data']) == 1:
				if data[key]['datatype'] == 'float':
					data[key]['data'] = float(data[key]['data'][0])
				else:
					data[key]['data'] = int(data[key]['data'][0])
	return data


def write_sorted(data):
	if not 'CellID' in data:
		exit('No CellID in ' + data['filename'])
	if data['CellID']['datatype'] != 'uint':
		exit('Unsupported datatype in ' + data['filename'])
	if data['CellID']['vectorsize'] != '1':
		exit('Unsupported vectorsize in ' + data['filename'])
	if data['CellID']['datasize'] != '4':
		exit('Unsupported datasize in ' + data['filename'])
	orig_cells = list(data['CellID']['data'])
	sort_cells = sorted(orig_cells)

	# data should be moved from old to new index
	mapping = dict()
	for old_i in range(len(orig_cells)):
		new_i = sort_cells.index(orig_cells[old_i])
		if old_i == new_i:
			continue
		mapping[old_i] = new_i

	# shuffle cell data around based on mapping
	variables = [key for key in data.keys() if 'type' in data[key] and data[key]['type'] == 'celldata']
	for var in variables:
		new_data = data[var]['data'].copy('K')
		for old_i in mapping.keys():
			new_data[mapping[old_i]] = data[var]['data'][old_i]
		data[var]['data'] = new_data

	# overwrite old data with new
	with open(data['filename'], 'r+b') as outfile:
		for var in variables:
			outfile.seek(data[var]['fileoffset'])
			outbytes = data[var]['data'].tobytes()
			outfile.write(outbytes)


if __name__ == '__main__':
	from sys import argv
	for filename in argv[1:]:
		write_sorted(load_data(load_metadata(filename)))
