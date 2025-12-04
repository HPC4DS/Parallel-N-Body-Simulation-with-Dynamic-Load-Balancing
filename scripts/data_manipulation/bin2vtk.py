# import os
#
# directory = "output/local/simulations/linear_barnes_hut_tree"
# f = directory + "/particles_iteration_00000.bin"
# print("filesize =", os.path.getsize(f))
# print("expected multiple =", 64)
# print("filesize % 64 =", os.path.getsize(f) % 64)

#==============================================================================

import numpy as np
import struct
import vtk
from vtkmodules.util.numpy_support import numpy_to_vtk
import glob

SPACE_DIM = 3
PARTICLE_SIZE = (
        SPACE_DIM * 8 +  # position
        SPACE_DIM * 8 +  # velocity
        8 +              # mass (double)
        8 +              # ID
        8                # morton_code (uint64)
)


STRUCT = struct.Struct("6d d Q Q")
# Meaning:
# 6 doubles (pos+vel)
# 1 double (mass)
# 1 size_t (ID)
# 1 uint64 (morton code)
# total = 9 doubles (72 bytes)

def load_particles(filename):
    with open(filename, "rb") as f:
        data = f.read()

    header_size = 4  # skip this
    particle_size = 72  # bytes per particle

    data = data[header_size:]  # remove header

    n = len(data) // particle_size
    if len(data) % particle_size != 0:
        raise RuntimeError("Corrupted file: remaining bytes = {}".format(len(data) % particle_size))

    arr = np.frombuffer(data, dtype=[
        ("pos",   np.float64, (3,)),
        ("vel",   np.float64, (3,)),
        ("mass",  np.float64),
        ("id",    np.uint64),
        ("code",  np.uint64),
    ])

    return arr

def save_vtp(particles, outname):
    points = vtk.vtkPoints()
    points.SetData(numpy_to_vtk(particles["pos"], deep=True))

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)

    # add velocity
    vel = numpy_to_vtk(particles["vel"], deep=True)
    vel.SetName("velocity")
    polydata.GetPointData().AddArray(vel)

    # add mass
    mass = numpy_to_vtk(particles["mass"], deep=True)
    mass.SetName("mass")
    polydata.GetPointData().AddArray(mass)

    # add ID
    pid = numpy_to_vtk(particles["id"], deep=True)
    pid.SetName("particle_id")
    polydata.GetPointData().AddArray(pid)

    # add Morton code
    code = numpy_to_vtk(particles["code"], deep=True)
    code.SetName("morton_code")
    polydata.GetPointData().AddArray(code)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(outname)
    writer.SetInputData(polydata)
    writer.SetDataModeToBinary()
    writer.Write()

# Process all
directory = "output/local/simulations/linear_barnes_hut_tree"
for f in sorted(glob.glob(directory + "/particles_iteration_*.bin")):
    frame = int(f.split("_")[-1].split(".")[0])
    out = f"{directory}/frame_{frame:05d}.vtp"
    print("writing", out)
    p = load_particles(f)
    save_vtp(p, out)
