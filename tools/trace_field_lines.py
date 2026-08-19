# field line tracing
import sys
from types import SimpleNamespace
import numpy as np
import scipy.interpolate as sci_intp
import scipy.integrate as sci_int
import rhybridpy as rhb
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from matplotlib.lines import Line2D
plt.rcParams.update({'font.size': 16})

save_figures = True # save figures as png files

Rp = 3390e3
input_file = 'state00004000.vlsv'
var_name = 'cellB' #'v_H+sw_ave'
box_margin_cells = 0.6 # outer stopping boundary inset, in grid cells
Rstop = 1.01*Rp # inner stopping radius [m]
empty_fill_value = np.nan

# integrator parameters
integ = SimpleNamespace()
integ.method = 'RK45'
integ.max_dx = 100e3 # maximum step size [m]
integ.rtol = 1e-3 # relative tolerance
integ.atol = 1e-6 # absolute tolerance
integ.S_MAX = 20*Rp # maximum arc length [m]
integ.intpol_method = 'linear' # interpolation method for integrated variable

# seed points
Npoints = 20

# circle orbit
Rc = 1.3*Rp
theta = np.linspace(0,2*np.pi,Npoints,endpoint=False)
x0 = Rc*np.cos(theta)
y0 = Rc*np.sin(theta)
z0 = np.zeros_like(x0)

# upstream region
#z0 = np.linspace(-2*Rp,2*Rp,Npoints)
#x0 = 2.0*Rp*np.ones_like(z0)
#y0 = 0.0*np.zeros_like(z0)

# gather seed points in one variable
r0 = np.column_stack((x0,y0,z0))
# handle single seed point
if r0.ndim == 1:
    r0 = r0.reshape(1,-1)
Npoints = r0.shape[0]
Rstop2 = Rstop * Rstop

# open file and check variable
vr = rhb.vlsvfile.VlsvReader(input_file)
if vr.check_variable(var_name) == False:
    print('ERROR: variable not found (' + var_name + ' / ' + input_file + ')')
    sys.exit()
var_dim = vr.read_variable_vectorsize(var_name)
if var_dim != 3:
    print('ERROR: variable should be a 3D vector (' + str(var_dim) + ' / ' + var_name + ')')
    sys.exit()

# simulation box dimensions
[xmin,ymin,zmin,xmax,ymax,zmax] = vr.get_spatial_mesh_extent()
[mx,my,mz] = vr.get_spatial_mesh_size() # how many blocks per direction
[sx,sy,sz] = vr.get_spatial_block_size() # how many cells per block per direction
nx = mx*sx # number of cells along x
ny = my*sy # number of cells along y
nz = mz*sz # number of cells along z
dx = (xmax-xmin)/nx
dy = (ymax-ymin)/ny
dz = (zmax-zmin)/nz
if not (dx == dy == dz):
    print(f'ERROR: grid cells must be cubic (dx={dx}, dy={dy}, dz={dz})')
    sys.exit()
# cell center coordinates
x_centers = xmin + (np.arange(nx) + 0.5)*dx
y_centers = ymin + (np.arange(ny) + 0.5)*dx
z_centers = zmin + (np.arange(nz) + 0.5)*dx

# read and reshape vector field
V = vr.read_variable(var_name)
cellids = vr.read_variable('CellID')
ncells = nx*ny*nz
if V.shape != (ncells, 3) or cellids.size != ncells:
    raise ValueError(
        f'Mesh/data size mismatch: expected {ncells} vector cells, '
        f'got V.shape={V.shape} and {cellids.size} CellIDs'
    )
cellids_sorted = np.argsort(cellids)
Vx = V[cellids_sorted,0].reshape(nz,ny,nx)
Vy = V[cellids_sorted,1].reshape(nz,ny,nx)
Vz = V[cellids_sorted,2].reshape(nz,ny,nx)

# vector field interpolators
interpolator_x = sci_intp.RegularGridInterpolator(
    (z_centers, y_centers, x_centers),
    Vx,
    bounds_error=False,
    fill_value=empty_fill_value,
    method=integ.intpol_method,
)
interpolator_y = sci_intp.RegularGridInterpolator(
    (z_centers, y_centers, x_centers),
    Vy,
    bounds_error=False,
    fill_value=empty_fill_value,
    method=integ.intpol_method,
)
interpolator_z = sci_intp.RegularGridInterpolator(
    (z_centers, y_centers, x_centers),
    Vz,
    bounds_error=False,
    fill_value=empty_fill_value,
    method=integ.intpol_method,
)

# get a unit vector of field orientation at x,y,z
def field_direction(t, crds):
    x, y, z = crds
    Vx = float(interpolator_x((z, y, x)))
    Vy = float(interpolator_y((z, y, x)))
    Vz = float(interpolator_z((z, y, x)))
    Vtot = np.sqrt(Vx*Vx + Vy*Vy + Vz*Vz)
    if not np.isfinite(Vtot):
        raise ValueError(f'Field interpolation is undefined at {crds}')
    if Vtot == 0:
        raise ValueError(f'Field direction is undefined at zero field: {crds}')
    return np.array([Vx,Vy,Vz])/Vtot

# outer box boundaries
def event_outer_boundaries(t, crds):
    x, y, z = crds
    return min(
        x - (xmin + dx*box_margin_cells),
        (xmax - dx*box_margin_cells) - x,
        y - (ymin + dx*box_margin_cells),
        (ymax - dx*box_margin_cells) - y,
        z - (zmin + dx*box_margin_cells),
        (zmax - dx*box_margin_cells) - z,
    )
# inner boundary
def event_inner_boundary(t, crds):
    x, y, z = crds
    r2 = x * x + y * y + z * z
    return r2 - Rstop2

event_outer_boundaries.terminal = True
event_outer_boundaries.direction = -1
event_inner_boundary.terminal = True
event_inner_boundary.direction = -1

# do field line tracing
print('integrator parameters:')
for name, value in vars(integ).items():
    print(name, '=', value)
print('tracing...')
field_lines_all = []
connectivity = []
total_forward_steps = 0
total_backward_steps = 0
for ii, seed in enumerate(r0, start=1):
    print('field line ' + str(ii) + "/" + str(Npoints))
    if np.linalg.norm(seed) <= Rstop:
        print(f'  skipped: seed is at or inside Rstop ({Rstop/Rp:.3f} Rp)')
        connectivity.append(3)
        field_lines_all.append(None)
        continue
    if event_outer_boundaries(0, seed) <= 0:
        print('  skipped: seed is outside the tracing box (' + str(seed/Rp) + ')')
        connectivity.append(3)
        field_lines_all.append(None)
        continue
    seed_as_float = seed.astype(float)
    try:
        # Forward tracing
        sol_forward = sci_int.solve_ivp(
            field_direction,
            [0, integ.S_MAX],
            seed_as_float,
            method=integ.method,
            rtol=integ.rtol,
            atol=integ.atol,
            max_step=integ.max_dx,
            events=[event_inner_boundary, event_outer_boundaries],
            vectorized=False,
        )
        # Backward tracing
        sol_backward = sci_int.solve_ivp(
            lambda t, Y: -field_direction(t, Y),
            [0, integ.S_MAX],
            seed_as_float,
            method=integ.method,
            rtol=integ.rtol,
            atol=integ.atol,
            max_step=integ.max_dx,
            events=[event_inner_boundary, event_outer_boundaries],
            vectorized=False,
        )
    except ValueError as exc:
        print(f'  tracing failed: {exc}')
        connectivity.append(3)
        field_lines_all.append(None)
        continue

    if not sol_forward.success or not sol_backward.success:
        print('  tracing failed: '
              f'forward={sol_forward.message}; backward={sol_backward.message}')
        connectivity.append(3)
        field_lines_all.append(None)
        continue

    total_forward_steps += len(sol_forward.t)
    total_backward_steps += len(sol_backward.t)
    # check and record field line connection type
    hit_fwd = sol_forward.t_events[0].size > 0
    hit_bwd = sol_backward.t_events[0].size > 0
    n_hits = int(hit_fwd) + int(hit_bwd)
    if n_hits == 2:
        conn = 2  # hit Rstop in both directions
    elif n_hits == 1:
        conn = 1   # hit Rstop only forward or only backward
    else:
        conn = 0      # never reached Rstop
    connectivity.append(conn)
    fwd_traj = sol_forward.y.T            # (n_fwd, 3)
    bwd_traj = sol_backward.y.T[::-1, :]  # (n_bwd, 3), reversed
    field_line = np.vstack((bwd_traj, fwd_traj[1:, :]))  # avoid duplicate seed
    field_lines_all.append(field_line)
print(f'Tracing complete. Total forward steps: {total_forward_steps}, '
      f'total backward steps: {total_backward_steps}.')

# plotting

ax = plt.figure(figsize=(10,10)).add_subplot(projection='3d')

# planet
u,v = np.meshgrid(np.linspace(0,np.pi,50),np.linspace(0,2*np.pi,50))
xs = 1 * np.sin(u) * np.cos(v)
ys = 1 * np.sin(u) * np.sin(v)
zs = 1 * np.cos(u)
cs = np.ones(xs.shape)
cs[xs > 0] = 0
planet_color_norm = colors.Normalize(vmin=0,vmax=1)
ax.plot_surface(xs,ys,zs,cmap=cm.Greys,facecolors=cm.Greys(planet_color_norm(cs)),shade=True, alpha=0.3, edgecolor='k')

colour_map_field_lines= {
    0: 'green',
    1: 'red',
    2: 'blue',
    3: 'gray',
    }
colours_field_lines = [colour_map_field_lines[c] for c in connectivity]
# field lines
for field_line, conn in zip(field_lines_all, connectivity):
    if field_line is None:
        continue
    colour = colour_map_field_lines[conn]
    ax.plot(field_line[:, 0]/Rp,field_line[:, 1]/Rp,field_line[:, 2]/Rp,'-', linewidth=1, alpha=0.8, c=colour)

# seed points
ax.scatter(r0[:,0]/Rp,r0[:,1]/Rp,r0[:,2]/Rp,color=colours_field_lines, marker='o', s=40)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_xlim(xmin / Rp, xmax / Rp)
ax.set_ylim(ymin / Rp, ymax / Rp)
ax.set_zlim(zmin / Rp, zmax / Rp)
ax.set_box_aspect([1, 1, 1])
ax.set_aspect('equal')
ax.set_title('field lines of ' + var_name)

legend_elements = [
    Line2D([0], [0], color='green',  lw=2, label='Not connected'),
    Line2D([0], [0], color='red', lw=2, label='Singly connected'),
    Line2D([0], [0], color='blue',   lw=2, label='Doubly connected'),
    Line2D([0], [0], color='gray', lw=2, label='Invalid / not traced'),
]

ax.legend(handles=legend_elements, loc='best', title='Connectivity')
plt.tight_layout()

if save_figures == True:
    plt.savefig('trace_field_lines.png')
    plt.clf()
    plt.close()
else:
    plt.show(block=True)
