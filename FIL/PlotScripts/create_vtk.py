#!/usr/bin/env python3
"""
Convert ETK HDF5 data to VTK time series for ParaView
Creates proper VTK structured grid files with time information
"""

import h5py
import numpy as np
import os
import glob
from collections import defaultdict

# Configuration
plane = "xy"  # Change to "xz" or "yz" if desired
refinement_level = 5  # Which refinement level to use
#data_dir = "data_hdf5_2D"
data_dir = "/mnt/raarchive/miler/FIL_RUNS/TOV_HOT_LONG/output-0000/data_hdf5_2D"
output_dir = "vtk_timeseries"

# Metric variables to exclude
EXCLUDE_VARS = [
    'alp', 'betax', 'betay', 'betaz',
    'gxx', 'gxy', 'gxz', 'gyy', 'gyz', 'gzz',
    'gtxx', 'gtxy', 'gtxz', 'gtyy', 'gtyz', 'gtzz',
    'W_z4c',
    'Bvec[0]', 'Bvec[1]', 'Bvec[2]'
]

# Create output directory
os.makedirs(output_dir, exist_ok=True)

def extract_variable_name(filename):
    """Extract variable name from filename"""
    basename = os.path.basename(filename)
    # Clean up variable names for VTK compatibility
    var_name = basename.replace(f'.{plane}.h5', '')
    var_name = var_name.replace('[', '_').replace(']', '')
    var_name = var_name.replace('-', '_')
    return var_name

def get_available_iterations(filename):
    """Get all available iterations from a file"""
    iterations = []
    with h5py.File(filename, 'r') as f:
        for key in f.keys():
            if 'it=' in key:
                try:
                    it = int(key.split('it=')[1].split()[0])
                    if it not in iterations:
                        iterations.append(it)
                except:
                    pass
    return sorted(iterations)

def load_data_for_iteration(filename, iteration, refinement_level):
    """Load and stitch data for a specific iteration"""
    
    with h5py.File(filename, 'r') as f:
        all_keys = list(f.keys())
        
        # Find components for this iteration and refinement level
        pattern_keys = [k for k in all_keys if f'it={iteration}' in k and f'rl={refinement_level}' in k]
        
        if not pattern_keys:
            return None, None, None, None
        
        # Get the full variable name
        full_var_name = pattern_keys[0].split(' it=')[0]
        pattern = f"{full_var_name} it={iteration} tl=0 rl={refinement_level}"
        components = [k for k in all_keys if k.startswith(pattern)]
        
        if not components:
            return None, None, None, None
        
        # Collect patches
        patches = []
        for comp in components:
            dataset = f[comp]
            data = dataset[:]
            attrs = dict(dataset.attrs)
            
            origin = attrs.get('origin', [0, 0, 0])
            delta = attrs.get('delta', [1, 1, 1])
            
            ny, nx = data.shape
            x = origin[0] + np.arange(nx) * delta[0]
            y = origin[1] + np.arange(ny) * delta[1]
            
            patches.append({
                'data': data,
                'x': x,
                'y': y
            })
        
        # Find global extent
        all_x = np.concatenate([p['x'] for p in patches])
        all_y = np.concatenate([p['y'] for p in patches])
        
        x_min, x_max = all_x.min(), all_x.max()
        y_min, y_max = all_y.min(), all_y.max()
        
        # Create uniform grid
        nx_grid = 400
        ny_grid = 400
        
        x_uniform = np.linspace(x_min, x_max, nx_grid)
        y_uniform = np.linspace(y_min, y_max, ny_grid)
        
        # Initialize with NaN
        data_grid = np.full((ny_grid, nx_grid), np.nan)
        
        # Fill from patches
        for patch in patches:
            x_patch = patch['x']
            y_patch = patch['y']
            data_patch = patch['data']
            
            ny_patch, nx_patch = data_patch.shape
            
            for iy in range(ny_patch):
                for ix in range(nx_patch):
                    ix_grid = np.argmin(np.abs(x_uniform - x_patch[ix]))
                    iy_grid = np.argmin(np.abs(y_uniform - y_patch[iy]))
                    data_grid[iy_grid, ix_grid] = data_patch[iy, ix]
        
        # Replace NaN with nearest neighbor
        mask = np.isnan(data_grid)
        if mask.any():
            from scipy.ndimage import distance_transform_edt
            indices = distance_transform_edt(mask, return_distances=False, return_indices=True)
            data_grid = data_grid[tuple(indices)]
        
        return x_uniform, y_uniform, data_grid, (x_min, x_max, y_min, y_max)

def write_vtk_structured_grid_binary(filename, x, y, z, scalar_data, scalar_names):
    """Write VTK structured grid file in BINARY format (more robust)"""
    
    nx = len(x)
    ny = len(y)
    nz = len(z)
    
    with open(filename, 'wb') as f:
        # Write ASCII header
        f.write(b"# vtk DataFile Version 3.0\n")
        f.write(b"ETK Data\n")
        f.write(b"BINARY\n")
        f.write(b"DATASET STRUCTURED_GRID\n")
        f.write(f"DIMENSIONS {nx} {ny} {nz}\n".encode('ascii'))
        f.write(f"POINTS {nx*ny*nz} float\n".encode('ascii'))
        
        # Write points in binary (big-endian floats)
        points = np.zeros((nx*ny*nz, 3), dtype='>f4')
        idx = 0
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    points[idx] = [x[i], y[j], z[k]]
                    idx += 1
        f.write(points.tobytes())
        
        # Write scalar data
        f.write(f"\nPOINT_DATA {nx*ny*nz}\n".encode('ascii'))
        
        for scalar_name, data in zip(scalar_names, scalar_data):
            f.write(f"\nSCALARS {scalar_name} float 1\n".encode('ascii'))
            f.write(b"LOOKUP_TABLE default\n")
            
            # Flatten and convert to big-endian float32
            scalar_flat = np.zeros(nx*ny*nz, dtype='>f4')
            idx = 0
            for k in range(nz):
                for j in range(ny):
                    for i in range(nx):
                        scalar_flat[idx] = data[j, i]
                        idx += 1
            
            f.write(scalar_flat.tobytes())

def write_pvd_file(vtk_files, iterations, output_file):
    """Write PVD file for time series"""
    
    with open(output_file, 'w') as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n')
        f.write('  <Collection>\n')
        
        for vtk_file, iteration in zip(vtk_files, iterations):
            basename = os.path.basename(vtk_file)
            f.write(f'    <DataSet timestep="{iteration}" file="{basename}"/>\n')
        
        f.write('  </Collection>\n')
        f.write('</VTKFile>\n')

# Main conversion
print(f"Converting ETK HDF5 to VTK time series")
print(f"Plane: {plane}, Refinement level: {refinement_level}")
print("=" * 60)

# Find all variables
pattern = os.path.join(data_dir, f"*.{plane}.h5")
files = sorted(glob.glob(pattern))

# Filter out excluded variables
variable_files = {}
for filename in files:
    var_name = extract_variable_name(filename)
    base_var = var_name.replace('_0', '').replace('_1', '').replace('_2', '')
    if base_var not in EXCLUDE_VARS:
        variable_files[var_name] = filename

print(f"Found {len(variable_files)} variables to convert (excluding metrics)")

# Get iterations from first file
first_file = list(variable_files.values())[0]
iterations = get_available_iterations(first_file)
print(f"Found {len(iterations)} iterations: {iterations}")

# Load all data
print("\nLoading data...")
all_data = defaultdict(dict)
x_coords = None
y_coords = None

for var_name, filename in variable_files.items():
    print(f"  Loading {var_name}...", end='', flush=True)
    
    for iteration in iterations:
        x, y, data, extent = load_data_for_iteration(filename, iteration, refinement_level)
        
        if data is not None:
            all_data[iteration][var_name] = data
            if x_coords is None:
                x_coords = x
                y_coords = y
    
    print(" ✓")

# Write VTK files
print("\nWriting VTK files...")
vtk_files = []
z_coords = np.array([0.0])  # 2D slice

for iteration in iterations:
    if iteration not in all_data or len(all_data[iteration]) == 0:
        continue
    
    output_file = os.path.join(output_dir, f"data_it{iteration:06d}.vtk")
    
    scalar_names = list(all_data[iteration].keys())
    scalar_data = [all_data[iteration][name] for name in scalar_names]
    
    write_vtk_structured_grid_binary(output_file, x_coords, y_coords, z_coords,
                                     scalar_data, scalar_names)
    
    vtk_files.append(output_file)
    print(f"  Written: {os.path.basename(output_file)} ({len(scalar_names)} variables)")

# Write PVD time series file
pvd_file = os.path.join(output_dir, "timeseries.pvd")
write_pvd_file(vtk_files, iterations, pvd_file)

print("\n" + "=" * 60)
print(f"Conversion complete!")
print(f"  Iterations: {len(vtk_files)}")
print(f"  Variables per iteration: {len(scalar_names)}")
print(f"\nTime series file: {pvd_file}")
print("\nTo view in ParaView:")
print(f"  1. Open ParaView")
print(f"  2. File → Open → {pvd_file}")
print(f"  3. Click 'Apply'")
print(f"  4. Use the time controls to animate through iterations")
