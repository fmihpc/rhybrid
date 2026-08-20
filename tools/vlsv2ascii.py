#! /usr/bin/env python3
'''
Converts RHybrid VLSV files to ASCII format.

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
		if name == 'VLSV':
			return
		if 'name' in attrs:
			name = attrs['name']
		current_name = name
		sim_data[name] = attrs

	def element_end(name):
		global current_name
		if name == 'VLSV':
			return
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
			if key == 'filename':
				continue

			if 'arraysize' in data[key] and data[key]['arraysize'] == 0:
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
					exit('only uint size 4 or 8 supported')
			elif data[key]['datatype'] == 'int':
				if data[key]['datasize'] == '8':
					dtype = 'int64'
				elif data[key]['datasize'] == '4':
					dtype = 'intc'
				else:
					exit('only int size 4 or 8 supported')
			if dtype == None:
				exit('unsupported datatype: ' + data[key]['datatype'])
			if data[key]['vectorsize'] != '1':
				dtype = data[key]['vectorsize'] + dtype

			infile.seek(data[key]['fileoffset'])
			data[key]['data'] = fromfile(
				infile, dtype = dtype, count = int(data[key]['arraysize']))

			if data[key]['arraysize'] == '1' and data[key]['vectorsize'] == '1':
				if data[key]['datatype'] == 'float':
					data[key]['data'] = float(data[key]['data'][0])
				else:
					data[key]['data'] = int(data[key]['data'][0])
	return data


# Saves cell data from given dict to ASCII file with _cells.txt suffix
def save_cell_data(data):
	cell_vars = sorted([
		key for key in data.keys() if
			'type' in data[key]
			and data[key]['type'] == 'celldata'
			and data[key]['arraysize'] != '1'
			and not key.startswith('MESH_NODE_CRDS_')])
	if len(cell_vars) == 0:
		return

	with open(data['filename'][:-5] + '_cells.txt', 'w') as outfile:
		outfile.write('# ASCII version of cell data in vlsv file '
			+ data['filename'] + '\n')
		for key in data.keys():
			if key == 'filename':
				continue
			if data[key]['arraysize'] != '1' or data[key]['vectorsize'] != '1':
				continue
			outfile.write('# {}: {}\n'.format(key, str(data[key]['data'])))

		nx = len(data['MESH_NODE_CRDS_X']['data']) - 1
		ny = len(data['MESH_NODE_CRDS_Y']['data']) - 1
		nz = len(data['MESH_NODE_CRDS_Z']['data']) - 1

		outfile.write('# Each cell on separate line:\n')
		outfile.write('# Columns 1,2,3: x,y,z components of cell center\n')
		col = 4
		for i in range(len(cell_vars)):
			var = cell_vars[i]
			if data[var]['vectorsize'] == '1':
				outfile.write('# Column {:d}: {}\n'.format(col, var))
				col += 1
			elif data[var]['vectorsize'] == '3':
				outfile.write('# Columns {:d},{:d},{:d}: x,y,z components of {}\n'
					.format(col, col+1, col+2, var))
				col += 3
			else:
				exit('Unsupported vectorsize: ' + data[var]['vectorsize'])

		# write cell data
		for z in range(nz):
			for y in range(ny):
				for x in range(nx):
					data_i = x + y*nx + z*nx*ny
					cell = data['CellID']['data'][data_i]
					cell_xi = cell % nx
					rx = 0.5 * (data['MESH_NODE_CRDS_X']['data'][cell_xi+1]
						+ data['MESH_NODE_CRDS_X']['data'][cell_xi])
					cell_yi = int(cell / nx) % ny
					ry = 0.5 * (data['MESH_NODE_CRDS_Y']['data'][cell_yi+1]
						+ data['MESH_NODE_CRDS_Y']['data'][cell_yi])
					cell_zi = int(cell / (nx * ny))
					rz = 0.5 * (data['MESH_NODE_CRDS_Z']['data'][cell_zi+1]
						+ data['MESH_NODE_CRDS_Z']['data'][cell_zi])
					outfile.write('{:.3e} {:.3e} {:.3e}'.format(rx, ry, rz))
					for var in cell_vars:
						fmt = ''
						if data[var]['datatype'] == 'float':
							fmt = ' {:.3e}'
						else:
							fmt = ' {:d}'
						if data[var]['vectorsize'] == '1':
							outfile.write(fmt.format(data[var]['data'][data_i]))
						elif data[var]['vectorsize'] == '3':
							dat = data[var]['data'][data_i]
							outfile.write((3*fmt).format(dat[0], dat[1], dat[2]))
					outfile.write('\n')

# Saves point data from given dict to ASCII file with _points.txt suffix
def save_point_data(data):
	point_vars = [key for key in data.keys() if
		'type' in data[key] and data[key]['type'] == 'point']
	if len(point_vars) > 1:
		exit('Unexpected number of point variables: ' + str(point_vars))
	if len(point_vars) == 0:
		return

	point_vars = point_vars[0]
	with open(data['filename'][:-5] + '_points.txt', 'w') as outfile:
		outfile.write('# ASCII version of point data in vlsv file '
			+ data['filename'] + '\n')
		for key in data.keys():
			if key == 'filename':
				continue
			if data[key]['arraysize'] != '1' or data[key]['vectorsize'] != '1':
				continue
			outfile.write('# {}: {}\n'.format(key, str(data[key]['data'])))

		if len(data[point_vars]['data']) == 0:
			outfile.write('# No point data in file ' + data['filename'] + '\n')
			return

		outfile.write('# Each point on separate line (lines ordered randomly):\n')
		outfile.write('# Columns 1,2,3: x,y,z components of point position\n')
		variables = sorted([
			key for key in data.keys() if
				'type' in data[key]
				and data[key]['type'] == 'pointdata'
				and not key.startswith('MESH_NODE_CRDS_')])
		col = 4
		for i in range(len(variables)):
			var = variables[i]
			if data[var]['vectorsize'] == '1':
				outfile.write('# Column {:d}: {}\n'.format(col, var))
				col += 1
			elif data[var]['vectorsize'] == '3':
				outfile.write('# Columns {:d},{:d},{:d}: x,y,z components of {}\n'
					.format(col, col+1, col+2, var))
				col += 3
			else:
				exit('Unsupported vectorsize: ' + data[var]['vectorsize'])

		# write point data
		for data_i in range(len(data[point_vars]['data'])):
			rx, ry, rz = data[point_vars]['data'][data_i]
			outfile.write('{:.5e} {:.5e} {:.5e}'.format(rx, ry, rz))
			for var in variables:
				fmt = ''
				if data[var]['datatype'] == 'float':
					fmt = ' {:.3e}'
				else:
					fmt = ' {:d}'
				if data[var]['vectorsize'] == '1':
					outfile.write(fmt.format(data[var]['data'][data_i]))
				elif data[var]['vectorsize'] == '3':
					dat = data[var]['data'][data_i]
					outfile.write((3*fmt).format(dat[0], dat[1], dat[2]))
			outfile.write('\n')


if __name__ == '__main__':
	from sys import argv
	for filename in argv[1:]:
		data = load_data(load_metadata(filename))
		save_cell_data(data)
		save_point_data(data)
