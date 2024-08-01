from dolfin import *
import numpy as np

class Io:

    def write_paraview(self, geometry, settings, file_name, fields, stamp=0):
        xdmf_file = XDMFFile(geometry.mesh.mpi_comm(), settings.result_folder + '/' + file_name + '.xdmf')

        xdmf_file.parameters["rewrite_function_mesh"] = False
        xdmf_file.parameters["functions_share_mesh"]  = True

        for field_name in fields:
            fields[field_name].rename(str(field_name), str(field_name))
            xdmf_file.write(fields[field_name], stamp)

        xdmf_file.close()