####################################################################################################
# Copyright (c) 2016 - 2022
# Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
#
# Author(s)
#       Marwan Abdellah <marwan.abdellah@epfl.ch>
#
# For complete list of authors, please see AUTHORS.md
#
# This file is part of Ultraliser <https://github.com/BlueBrain/Ultraliser>
#
# This library is free software; you can redistribute it and/or modify it under the terms of the
# GNU Lesser General Public License version 3.0 as published by the Free Software Foundation.
#
# This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# See the GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along with this library;
# if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301 USA.
# You can also find it on the GNU web site < https://www.gnu.org/licenses/gpl-3.0.en.html >
####################################################################################################


# System imports
import argparse
import nrrd
import struct
from pathlib import Path

from pandas import wide_to_long

####################################################################################################
# @parse_command_line_arguments
####################################################################################################
def parse_command_line_arguments(arguments=None):

    # Add all the options
    description = 'Creating Ultraliser-compatible RAW/HDR volumes from NRRD files using pynrrd'
    parser = argparse.ArgumentParser(description=description)

    arg_help = 'Input NRRD file'
    parser.add_argument('--input-file', action='store', dest='input_file', help=arg_help)

    arg_help = 'Output directory where output written'
    parser.add_argument('--output-directory', action='store', dest='output_directory', help=arg_help)
                        
    # Parse the arguments
    return parser.parse_args()

####################################################################################################
# @write_raw_header_file
####################################################################################################
def write_raw_header_file(prefix, file_type, width, height, depth):

    # Open the header file 
    fhandle = open('%s.hdr' % prefix, 'w')

    # File type 
    fhandle.write('%s \n' % file_type)

    # Dimensions 
    fhandle.write('%d ' % width)
    fhandle.write('%d ' % height)
    fhandle.write('%d\n' % depth)

    # Close the file 
    fhandle.close()


####################################################################################################
# @write_raw_file
####################################################################################################
def write_raw_file(prefix, data, width, height, depth):

    # Open the data file 
    fhandle = open('%s.img' % prefix, 'wb')

    # Write the data
    for i in range(0, width):
        print(i, width)
        for j in range(0, height):
            for k in range(0, depth):
                fhandle.write(struct.pack('>H', data[i][j][k]))
                
    
    # Close the file 
    fhandle.close()

def write_nrrd_file(prefix, data, width, height, depth):

     # Open the data file 
    fhandle = open('%s.nrrd' % prefix, 'w')
    fhandle.write('dimension: 3\n')
    fhandle.write('type: uint16\n')
    fhandle.write('sizes: %s %s %s\n' % (str(width), str(height), str(depth)))
    fhandle.write('endian: little\n')
    fhandle.write('encoding: raw\n')
    fhandle.close()

    fhandle = open('%s.nrrd' % prefix, 'ab')
    
    # Write the data
    for i in range(0, width):
        print(i, width - 1)
        for j in range(0, height):
            for k in range(0, depth):
                buffer = struct.pack('H', data[i][j][k])
                fhandle.write(buffer)

    # Close the file 
    fhandle.close()


####################################################################################################
# @parse_command_line_arguments
####################################################################################################
def convert_nrrd_to_raw(arguments):

    # Get the file name from the 
    file_name = Path(arguments.input_file).stem

    # Construct the prefix 
    prefix = '%s/%s' % (arguments.output_directory, file_name)

    # Read the NRRD file 
    nrrd_data, nrrd_header = nrrd.read(arguments.input_file)
    
    # Get the NRRD type 
    nrrd_type = nrrd_header['type']

    # Get the NRRD dimensions 
    nrrd_dimension = nrrd_header['dimension']

    if (nrrd_dimension != 3):
        print('NRRD dimensions are less than 3, cannot be a volume! EXITTING')
        exit(0)

    # Get the volume dimensions 
    volume_dimensions = nrrd_header['sizes']
    volume_width = volume_dimensions[0]
    volume_height = volume_dimensions[1]
    volume_depth = volume_dimensions[2]

    if (nrrd_type == "signed char" or 
        nrrd_type == "int8" or 
        nrrd_type ==  "int8_t"):

        # Not supported 
        print('File type NOT supported: EXITTING!')
        exit(0)
        

    elif (nrrd_type == "uchar" or 
          nrrd_type == "unsigned char" or 
          nrrd_type == "uint8" or 
          nrrd_type == "uint8_t"):

        write_raw_header_file(prefix=prefix, file_type='u8', 
            width=volume_width, height=volume_height, depth=volume_depth)
        
        write_raw_file(prefix=prefix, data=nrrd_data, 
            width=volume_width, height=volume_height, depth=volume_depth)


    elif (nrrd_type == "short" or 
          nrrd_type == "short int" or 
          nrrd_type == "signed short" or 
          nrrd_type == "signed short int" or 
          nrrd_type == "int16" or 
          nrrd_type == "int16_t"):
        
        # Not supported 
        print('File type NOT supported: EXITTING!')
        exit(0)

    elif (nrrd_type == "ushort" or 
          nrrd_type == "short int" or 
          nrrd_type == "unsigned short" or 
          nrrd_type == "unsigned short int" or 
          nrrd_type == "uint16" or 
          nrrd_type == "uint16_t"):
        
        write_raw_header_file(prefix=prefix, file_type='u16', 
            width=volume_width, height=volume_height, depth=volume_depth)

        #write_raw_file(prefix=prefix, data=nrrd_data, 
        #    width=volume_width, height=volume_height, depth=volume_depth)

        write_nrrd_file(prefix=prefix, data=nrrd_data, 
            width=volume_width, height=volume_height, depth=volume_depth)


    elif (nrrd_type == "int" or  
          nrrd_type == "signed int" or
          nrrd_type == "int32" or 
          nrrd_type == "int32_t"):

        # Not supported 
        print('File type NOT supported: EXITTING!')
        exit(0)
        
    elif (nrrd_type == "uint" or  
          nrrd_type == "unsigned int" or
          nrrd_type == "uint32" or 
          nrrd_type == "uint32_t"):

        write_raw_header_file(prefix=prefix, file_type='u32', 
            width=volume_width, height=volume_height, depth=volume_depth)

        write_raw_file(prefix=prefix, data=nrrd_data, 
            width=volume_width, height=volume_height, depth=volume_depth)


    elif (nrrd_type == "longlong" or 
          nrrd_type == "long long" or 
          nrrd_type == "long long int" or 
          nrrd_type == "signed long long" or 
          nrrd_type == "signed long long int" or 
          nrrd_type == "int64" or 
          nrrd_type == "int64_t"):

        # Not supported 
        print('File type NOT supported: EXITTING!')
        exit(0)

    elif (nrrd_type == "ulonglong" or 
          nrrd_type == "unsigned long long" or 
          nrrd_type == "unsigned long long int" or 
          nrrd_type == "uint64" or 
          nrrd_type == "uint64_t"):
        
        write_raw_header_file(prefix=prefix, file_type='u64', 
            width=volume_width, height=volume_height, depth=volume_depth)

        write_raw_file(prefix=prefix, data=nrrd_data, 
            width=volume_width, height=volume_height, depth=volume_depth)

    
    elif (nrrd_type == "float"):
        
        write_raw_header_file(prefix=prefix, file_type='f32', 
            width=volume_width, height=volume_height, depth=volume_depth) 

        write_raw_file(prefix=prefix, data=nrrd_data, 
            width=volume_width, height=volume_height, depth=volume_depth)

    elif (nrrd_type == "double"):
        
        write_raw_header_file(prefix=prefix, file_type='f64', 
            width=volume_width, height=volume_height, depth=volume_depth) 

        write_raw_file(prefix=prefix, data=nrrd_data, 
            width=volume_width, height=volume_height, depth=volume_depth)

        
    elif (nrrd_type == "block"):
        
        # Not supported 
        print('File type NOT supported: EXITTING!')
        exit(0)

    else:
        print("Cannot recognize the file type! EXITTING")
        exit(0) 

####################################################################################################
# @ Main
####################################################################################################
if __name__ == "__main__":

    # Parse the command line arguments
    args = parse_command_line_arguments()

    # Do it
    convert_nrrd_to_raw(arguments=args)
