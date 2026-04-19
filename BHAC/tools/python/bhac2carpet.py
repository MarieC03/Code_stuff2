import sys
import numpy as np
import h5py
import multiprocessing as mp

###### Functions ######################################

def make_uniform(dset, Levels, LogLocs, N_blocks, extentX1, extentX2, extentX3, GZ, max_ref, nxlone):
    max_ref -= 1
    nx1 = extentX1 - 2*GZ
    nx2 = extentX2 - 2*GZ
    nx3 = extentX3 - 2*GZ

    Nx1 = nx1 * 2**max_ref * nxlone[0] / (extentX1-2*GZ)
    Nx2 = nx2 * 2**max_ref * nxlone[1] / (extentX2-2*GZ)
    Nx3 = nx3 * 2**max_ref * nxlone[2] / (extentX3-2*GZ)
    new = np.zeros((Nx3, Nx2, Nx1))

    LogLocsX = LogLocs[:,0] - 1
    LogLocsY = LogLocs[:,1] - 1
    LogLocsZ = LogLocs[:,2] - 1

    for i in range(0,N_blocks):
        mylevel = Levels[i]
        refs = max_ref - mylevel + 1
        tmp = dset[i,GZ:-GZ,GZ:-GZ,GZ:-GZ]
        for c in range(0,tmp[0,0,:].size):
             for b in range(0,tmp[0,:,0].size):
                 for a in range(0,tmp[:,0,0].size):
                     amin = a * 2**refs + LogLocsX[i]*(extentX1-2*GZ) * 2**refs
                     amax = amin + 2**refs
                     bmin = b * 2**refs + LogLocsY[i]*(extentX2-2*GZ) * 2**refs
                     bmax = bmin + 2**refs
                     cmin = c * 2**refs + LogLocsZ[i]*(extentX3-2*GZ) * 2**refs
                     cmax = cmin + 2**refs

                     new[cmin:cmax,bmin:bmax,amin:amax] = tmp[c,b,a]

    return new

##############################################################

def find_dx(xs,Levels,max_ref):
    for ref in range(1,max_ref):
        for i,L in enumerate(Levels):
            if(L == ref):
                return abs(xs[i,1]-xs[i,0])*0.5**(max_ref - ref)

#############################################################

def get_nxlone(xs,GZ,Levels,max_ref):
    xmaxGZ = np.max(xs)
    xminGZ = np.min(xs)
    dxlone = find_dx(xs, Levels, max_ref) / 0.5**(max_ref - 1)
    xmax = xmaxGZ - GZ*dxlone
    xmin = xminGZ + GZ*dxlone

    return int((xmax - xmin) / dxlone + 1)

#############################################################

def compute_dset(index):

    filename = bhac_basename + str(index).zfill(4) + ".hdf5"
    file = h5py.File(filename, "r")

    dset = file["Cell-centered Variables"]
    dset = dset[:,int(v_index),:,:,:]

    Levels = file["Levels"]
    LogLocs = file["LogicalLocations"]

    N_blocks = file.attrs["Nleafs"]
    it = file.attrs["Iteration"]
    time = file.attrs["Time"]

    dset = make_uniform(dset, Levels, LogLocs, N_blocks, extentX1, extentX2, extentX3, GZ, max_ref, nxlone)

    return (it, time, dset)

#############################################################


if ((sys.argv[1] == "-h") or (sys.argv[1] == "--help")):
    print """This script converts BhacHDF5 data to hdf5 data in the Carpet format. This 
is especially useful for 3D data visualisation with Amira. It is therefore assumed
that the data to be converted is 3D. 1D and 2D data conversion is not yet supported.

For AMR data the script will produce a Carpet file with only one refinement level,
which is going to be the highest level from the input data. For coordinate systems
other than Cartesian data will be interpolated onto a Cartesian grid as this is the
only format supported by Carpet and hence AMIRA (this interpolation is not yet
implemented.)

Usage: python2 bhac2carpet.py <InputBhacHDF5_Basename> <bhac_index_start> <bhac_index_end> <variable_index> <OutputCarpetHDF5File.h5>"""
    sys.exit()

elif (len(sys.argv) != 6):
    print "Usage: python2 bhac2carpet.py <InputBhacHDF5_Basename> <bhac_index_start> <bhac_index_end> <variable_index> <OutputCarpetHDF5File.h5>"
    sys.exit()

bhac_basename = sys.argv[1]
if (bhac_basename[-1:-5:-1] == "5fdh"):
    print "Only type BhacHDF5_Basename, e.g. output/data0003.hdf5 --> basename: output/data"
    sys.exit()

start = sys.argv[2]
end = sys.argv[3]
v_index = sys.argv[4]
carpet_filename = sys.argv[5]

print "Reading variable ", v_index, " from ", bhac_basename, start, "-", end

carpet = h5py.File(carpet_filename, "w")

g = carpet.create_group("Parameters and Global Attributes")
g.attrs.create("nioprocs", 1, dtype="int32")

#Some parameters that are the same for all datasets
name = "BHAC::data"

filename = bhac_basename + start + ".hdf5"
file = h5py.File(filename, "r")

extentX1 = file.attrs["ExtentX1"]
extentX2 = file.attrs["ExtentX2"]
extentX3 = file.attrs["ExtentX3"]
GZ = file.attrs["GhostZones"]
max_ref = file.attrs["LevMax"]

Levels = file["Levels"]
xs1 = file["x1v"]
xs2 = file["x2v"]
xs3 = file["x3v"]

dx = find_dx(xs1, Levels, max_ref)
dy = find_dx(xs2, Levels, max_ref)
dz = find_dx(xs3, Levels, max_ref)

originX = xs1[0,GZ]
originY = xs2[0,GZ]
originZ = xs3[0,GZ]

nxlone = [0,0,0]
nxlone[0] = get_nxlone(xs1,GZ,Levels,max_ref)
nxlone[1] = get_nxlone(xs2,GZ,Levels,max_ref)
nxlone[2] = get_nxlone(xs3,GZ,Levels,max_ref)

pool = mp.Pool(processes=4)
output = [pool.apply_async(compute_dset, args=(i,)) for i in range(int(start),int(end))]
results = [p.get() for p in output]

print "Creating Carpet file, carpet_filename"
#Create Carpet dataset from results
for i in range(int(start),int(end)):
    it = results[i][0]
    time = results[i][1]
    dset = results[i][2]
    carpet_dset = carpet.create_dataset("%s it=%d tl=0 m=0 rl=0 c=0" % (name, it),
                         shape=dset.shape, dtype=dset.dtype, data=dset) 
    
    carpet_dset.attrs["origin"] = np.array([originX, originY, originZ], dtype="float")
    carpet_dset.attrs["iorigin"] = np.array([0]*3,dtype="int32")
    carpet_dset.attrs.create("level", 0, dtype="int32")
    carpet_dset.attrs.create("timestep", it, dtype="int32")
    carpet_dset.attrs.create("time", time, dtype="float")
    carpet_dset.attrs["delta"] = np.array([dx,dy,dz], dtype="float")
    carpet_dset.attrs["name"] = np.string_(name + "\0")

print "All done."

