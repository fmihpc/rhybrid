# 
# This file is part of Analysator.
# Copyright 2013-2025 Finnish Meteorological Institute
# Copyright 2017-2018 University of Helsinki
# 
# For details of usage, see the COPYING file and read the "Rules of the Road"
# at http://www.physics.helsinki.fi/vlasiator/
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# 

import logging
import struct
import xml.etree.ElementTree as ET
import ast
import numpy as np
import os
from reduction import datareducers,data_operators
from vlsvwriter import VlsvWriter
from variable import get_data

def get_numpy_data_type(dataType,dataSize):
   ''' Get numpy data type corresponding to datatype and datasize of
       a VLSV file variable
   '''
   if dataType == "float" and dataSize == 8:
      return np.float64
   elif dataType == "float" and dataSize == 4:
      return np.float32
   elif dataType == "uint" and dataSize == 4:
      return np.uint32
   elif dataType == "uint" and dataSize == 8:
      return np.uint64

class VlsvParticles(object):
   ''' Class for reading point mesh based particle outputs from RHybrid
   '''
   file_name=""
   def __init__(self, file_name):
      ''' Initializes the vlsv file (opens the file, reads the file footer and reads in some parameters)
          :param file_name:     Name of the vlsv file
      '''
      # Make sure the path is set in file name: 
      file_name = os.path.abspath(file_name)

      self.file_name = file_name
      self.__fptr = open(self.file_name,"rb")
      self.__xml_root = ET.fromstring("<VLSV></VLSV>")
      self.__fileindex_for_cellid={}
      self.__fileindex_for_cellid_blocks={}
      self.__read_xml_footer()
      self.__fptr.close()

   def __read_xml_footer(self):
      ''' Reads in the XML footer of the VLSV file and store all the content
      '''
      max_xml_size = 1000000
      #(endianness,) = struct.unpack("c", fptr.read(1))
      if self.__fptr.closed:
         fptr = open(self.file_name,"rb")
      else:
         fptr = self.__fptr
      # Eight first bytes indicate whether the system is big_endianness or something else
      endianness_offset = 8
      fptr.seek(endianness_offset)
      # Read 8 bytes as unsigned long long (uint64_t in this case) after endianness, this tells the offset of the XML file.
      uint64_byte_amount = 8
      (offset,) = struct.unpack("Q", fptr.read(uint64_byte_amount))
      # Move to the xml offset
      fptr.seek(offset)
      # Read the xml data
      xml_data = fptr.read(max_xml_size)
      # Read the xml as string
      (xml_string,) = struct.unpack("%ds" % len(xml_data), xml_data)
      # Input the xml data into xml_root
      self.__xml_root = ET.fromstring(xml_string)
      if self.__fptr.closed:
         fptr.close()

   def list(self):
      ''' Print out the particle point mesh content of the file. Useful
          for interactive usage.
      '''
      logging.info("tag = MESH, type = point")
      for child in self.__xml_root:
         if child.tag == "MESH" and "type" in child.attrib:
             if child.attrib["type"] == "point":
                 logging.info("   " + str(child.attrib["name"]))
      logging.info("tag = VARIABLE, type = pointdata")
      for child in self.__xml_root:
         if child.tag == "VARIABLE" and "type" in child.attrib:
             if child.attrib["type"] == "pointdata":
                 logging.info("   " + str(child.attrib["name"]))

   def get_all_point_meshes(self):
      ''' Returns names of all particle point meshes in the file
          :returns: List of particle point mesh names in the file
          .. code-block:: python
             # Example usage:
             vrp = pt.vlsvparticle.VlsvParticle("test.vlsv")
             names_pmeshes = vrp.get_all_point_meshes()
      '''
      meshList = [];
      for child in self.__xml_root:
         if child.tag == "MESH" and "type" in child.attrib:
             if child.attrib["type"] == "point":
                 meshList.append(child.attrib["name"])
      return meshList

   def get_all_point_variables(self):
      ''' Returns names of all particle point mesh variables in the file
          :returns: List of particle point mesh variable names in the file
          .. code-block:: python
             # Example usage:
             vrp = pt.vlsvparticle.VlsvParticle("test.vlsv")
             names_pvars = vrp.get_all_point_variables()
      '''
      varList = [];
      for child in self.__xml_root:
         if child.tag == "VARIABLE" and "type" in child.attrib:
             if child.attrib["type"] == "pointdata":
                 varList.append(child.attrib["name"])
      return varList

   def read_particle_mesh(self,name):
      ''' Read a particle point mesh
          Arguments:
          :param name: Name of the mesh
          :returns: numpy array with the x, y and z coordinates of the particles
          .. code-block:: python
             # Example usage:
             vrp = pt.vlsvparticle.VlsvParticle("test.vlsv")
             ples = vrp.read_particle_mesh("ParticlePointMesh_sw_H+")
      '''
      if self.__fptr.closed:
         fptr = open(self.file_name,"rb")
      else:
         fptr = self.__fptr

      res = None
      for child in self.__xml_root:
          if child.tag == "MESH" and "type" in child.attrib:
              if child.attrib["type"] != "point" or child.attrib["name"] != name:
                  continue
              vector_size = ast.literal_eval(child.attrib["vectorsize"])
              if vector_size != 3:
                  logging.info(" WARNING, number of coordinates of a point mesh not three: " + name + ": " + str(vector_size))
              array_size = ast.literal_eval(child.attrib["arraysize"])
              element_size = ast.literal_eval(child.attrib["datasize"])
              datatype = child.attrib["datatype"]
              offset = ast.literal_eval(child.text)
              fptr.seek(offset)
              res = np.fromfile(fptr, dtype=self.get_numpy_data_type(datatype,element_size), count=vector_size*array_size)
              res = res.reshape(array_size, vector_size)
              return res
      if self.__fptr.closed:
         fptr.close()
      logging.info(" WARNING, particle point mesh not found: " + name)
      return res

   def read_particle_variable(self,name):
      ''' Read a particle point variable
          Arguments:
          :param name: Name of the variable
          :returns: numpy array with the data
          .. code-block:: python
             # Example usage:
             vrp = pt.vlsvparticle.VlsvParticle("test.vlsv")
             v = vrp.read_particle_variable("v_sw_H+")
      '''
      if self.__fptr.closed:
         fptr = open(self.file_name,"rb")
      else:
         fptr = self.__fptr

      res = None
      for child in self.__xml_root:
          if child.tag == "VARIABLE" and "type" in child.attrib:
              if child.attrib["type"] != "pointdata" or child.attrib["name"] != name:
                  continue
              vector_size = ast.literal_eval(child.attrib["vectorsize"])
              array_size = ast.literal_eval(child.attrib["arraysize"])
              element_size = ast.literal_eval(child.attrib["datasize"])
              datatype = child.attrib["datatype"]
              offset = ast.literal_eval(child.text)
              fptr.seek(offset)
              res = np.fromfile(fptr, dtype=get_numpy_data_type(datatype,element_size), count=vector_size*array_size)
              res = res.reshape(array_size, vector_size)
              return res
      if self.__fptr.closed:
         fptr.close()
      logging.info(" WARNING, particle variable not found: " + name)
      return res

