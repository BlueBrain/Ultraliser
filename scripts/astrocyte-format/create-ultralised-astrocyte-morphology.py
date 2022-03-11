# System imports
import argparse
import h5py
from pathlib import Path
from archngv import NGVCircuit
import pandas as pd



####################################################################################################
# CONSTANTS
####################################################################################################
# The directory that stores morphology points in an .H5 file
H5_POINTS_DIRECTORY = '/points'

# The directory that stores morphology connectivity in an .H5 file
H5_STRUCTURE_DIRECTORY = '/structure'

# The directory where the endfeet points will be stored
H5_ENDFEET_POINTS_DIRECTORY = '/endfeetpoints'

# The directory where the endfeet patches will be stored
H5_ENDFEET_PATCHES_DIRECTORY = '/endfeetpoints'

# The directory that stores morphology perimeters in an .H5 file
# This is only used for glia, but not used for neurons
H5_PERIMETERS_DIRECTORY = '/perimeters'


####################################################################################################
# @parse_command_line_arguments
####################################################################################################
def parse_command_line_arguments(arguments=None):

    # Add all the options
    description = 'Creating Ultraliser-compatible astrocyte morphology with endfeet'
    parser = argparse.ArgumentParser(description=description)

    arg_help = 'BBP Circuit'
    parser.add_argument('--circuit', action='store', dest='circuit', help=arg_help)

    arg_help = 'Astrocyte GID 0'
    parser.add_argument('--gid-0', action='store', dest='gid_0', type=int, help=arg_help)

    arg_help = 'Astrocyte GID N'
    parser.add_argument('--gid-n', action='store', dest='gid_n', type=int, help=arg_help)

    arg_help = 'Output directory where output written'
    parser.add_argument('--output-directory', action='store', dest='output_directory', help=arg_help)
                        
    # Parse the arguments
    return parser.parse_args()


####################################################################################################
# @EndFoot
####################################################################################################
class EndFoot:

    ################################################################################################
    # @__init__
    ################################################################################################
    def __init__(self, vertices=None, triangles=None, thickness=None):
        """Constructor
        """
        self.vertices = vertices
        self.triangles = triangles
        self.thickness = thickness


####################################################################################################
# @create_full_astrocyte_structure
####################################################################################################
def create_full_astrocyte_structure(astrocyte_h5_file,
                                    coordinates,
                                    endfeet,
                                    gid,
                                    output_directory):

    # Load h5 file and read its points and connectivity
    astrocyte_skeleton = h5py.File(astrocyte_h5_file, 'r')

    # Read the points
    points_list = astrocyte_skeleton[H5_POINTS_DIRECTORY]

    # Read the structure
    structure_list = astrocyte_skeleton[H5_STRUCTURE_DIRECTORY]

    # Create a new h5 file and write the directories
    full_astrocyte_h5_file = h5py.File('%s/%d.h5' % (output_directory, gid), "w")

    # Copy the points dataset
    full_astrocyte_h5_file.create_dataset(name="points",
                                          shape=points_list.shape,
                                          dtype=points_list.dtype,
                                          data=points_list)

    # Copy the structure list
    full_astrocyte_h5_file.create_dataset(name="structure",
                                          shape=structure_list.shape,
                                          dtype=structure_list.dtype,
                                          data=structure_list)

    # Add the coordinates
    full_astrocyte_h5_file.create_dataset(name="coordinates",
                                          shape=(1, 3),
                                          dtype='float64',
                                          data=coordinates)

    # Create the endfeet points (or vertices)
    vertex_start_index = 0
    vertex_last_index = -1
    vertex_indices = list()
    vertices_data = list()

    # Triangles
    triangle_start_index = 0
    triangle_last_index = -1
    triangle_indices = list()
    triangle_data = list()

    # Collect the vertices from the endfeet data
    for eid, endfoot in enumerate(endfeet):

        vertex_start_index = vertex_last_index + 1
        vertex_last_index = vertex_start_index + len(endfoot.vertices)
        vertex_indices.append([vertex_start_index, vertex_last_index])

        triangle_start_index = triangle_last_index + 1
        triangle_last_index = triangle_start_index + len(endfoot.triangles)
        triangle_indices.append([triangle_start_index, triangle_last_index])

        # Compose the lists
        for vid, vertex in enumerate(endfoot.vertices):
            vertices_data.append([vertex[0], vertex[1], vertex[2], endfoot.thickness])

        # Compose the lists
        for tid, triangle in enumerate(endfoot.triangles):
            triangle_data.append([triangle[0], triangle[1], triangle[2]])

    full_astrocyte_h5_file.create_dataset(name="endfeet_vertex_indices",
                                          shape=(len(vertex_indices), 2),
                                          dtype='int32',
                                          data=vertex_indices)

    full_astrocyte_h5_file.create_dataset(name="endfeet_vertex_data",
                                          shape=(len(vertices_data), 4),
                                          dtype='float64',
                                          data=vertices_data)

    full_astrocyte_h5_file.create_dataset(name="endfeet_triangle_indices",
                                          shape=(len(triangle_indices), 2),
                                          dtype='int32',
                                          data=triangle_indices)

    full_astrocyte_h5_file.create_dataset(name="endfeet_triangle_data",
                                          shape=(len(triangle_data), 3),
                                          dtype='int64',
                                          data=triangle_data)


####################################################################################################
# @create_ultralised_astrocytes
####################################################################################################
def create_ultralised_astrocytes(circuit_path,
                                 gid_0,
                                 gid_n,
                                 output_directory):
    # Load the circuit
    circuit = NGVCircuit(circuit_path)

    # Load astrocyte data
    astrocytes = circuit.astrocytes

    # Load the connectivity
    gv_conn = circuit.gliovascular_connectome

    # Get the positions and convert them into numpy arrays
    astrocyte_positions = astrocytes.positions().to_numpy()

    # Get the surface meshes
    surface_meshes = gv_conn.surface_meshes

    from joblib import Parallel, delayed
    Parallel(n_jobs=args.number_cores)(delayed(run_command)(i) for i in commands)


    # For the astrocytes in range
    for i in range(gid_0, gid_n + 1):

        # Get the IDs of the endfeet
        endfeet_ids = gv_conn.astrocyte_endfeet(i)

        # Get the skeleton H5 file
        astrocyte_skeleton_h5_file = astrocytes.morph.get_filepath(i)

        # Get the position of the astrocyte
        astrocyte_coordinates = list(astrocyte_positions[i])

        # A list of empty endfeet
        endfeet = list()
        for eid in endfeet_ids:
            endfoot = EndFoot(vertices=surface_meshes.mesh_points(eid),
                              triangles=surface_meshes.mesh_triangles(eid),
                              thickness=surface_meshes.get("surface_thickness", eid))

            # Append the endfoot to the list
            endfeet.append(endfoot)

        # Create the ultralised astrocyte
        create_full_astrocyte_structure(astrocyte_h5_file=astrocyte_skeleton_h5_file,
                                        coordinates=astrocyte_coordinates,
                                        endfeet=endfeet,
                                        gid=i,
                                        output_directory=output_directory)


####################################################################################################
# @ Main
####################################################################################################
if __name__ == "__main__":

    # Parse the command line arguments
    args = parse_command_line_arguments()

    # Do it
    create_ultralised_astrocytes(circuit_path=args.circuit,
                                 gid_0=args.gid_0,
                                 gid_n=args.gid_n,
                                 output_directory=args.output_directory)



