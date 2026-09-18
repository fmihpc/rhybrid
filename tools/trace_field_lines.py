"""Trace 3-D field lines on a uniform, cell-centred RHybrid grid.

The default backend uses a Numba-compiled trilinear interpolator and normalizes
the interpolated vector inside the same compiled function. The vector-valued
SciPy interpolator remains available as a slower reference backend.

Usage:
    python trace_field_lines.py

Connectivity values used for plotting are 0=not connected, 1=singly connected,
2=doubly connected, and 3=invalid, skipped, or length-limited.
"""

import os
import sys
from time import perf_counter
from types import SimpleNamespace
import numpy as np
import scipy.interpolate as sci_intp
import scipy.integrate as sci_int
import rhybridpy as rhb
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from matplotlib.lines import Line2D

numba_import_error = None
try:
    from numba import njit
except (ImportError, AttributeError) as exc:
    # AttributeError also covers an incompatible Numba/coverage installation.
    njit = None
    numba_import_error = exc

plt.rcParams.update({'font.size': 16})

def print_timing(label, start):
    """Print elapsed wall-clock time since *start* and return a new start."""
    now = perf_counter()
    print(f'TIMING {label}: {now - start:.3f} s')
    return now

# user configuration

Rp = 3390e3 # planet radius
file_input_state = './state00004000.vlsv';
file_output_figure = 'trace_field_lines.png'
save_output_figure = True
var_name = 'cellB' # vector variable to be traced
box_margin_cells = 1.0 # outer stopping boundary inset in dx, in grid cells (keep safely inside the outermost interpolation centres)
Rstop = 1.01*Rp # inner stopping radius [m]

# integrator parameters
integ = SimpleNamespace()
integ.method = 'RK45'
integ.max_step_cells = 0.1     # maximum RK45 step, measured in grid-cell widths
integ.error_cells = 1e-4       # approximate local position-error target, in grid-cell widths
integ.max_length = 20*Rp       # maximum arc length traced in each direction [m]
integ.intpol_method = 'linear' # interpolation method for integrated variable
integ.intpol_backend = 'scipy' # 'numba' (fast) or 'scipy' (reference)

if integ.intpol_backend not in ('numba', 'scipy'):
    raise ValueError(
        f"integ.intpol_backend must be 'numba' or 'scipy', "
        f'got {integ.intpol_backend!r}'
    )
if integ.intpol_backend == 'numba' and njit is None:
    raise RuntimeError(
        'Numba is required when integ.intpol_backend is set to numba. '
        f'Install/fix numba or select scipy. Import error: {numba_import_error}'
    )
if integ.intpol_backend == 'numba' and integ.intpol_method != 'linear':
    raise ValueError(
        "The Numba backend currently supports only linear interpolation"
    )
for parameter_name in ('max_step_cells', 'error_cells', 'max_length'):
    parameter_value = getattr(integ, parameter_name)
    if not np.isfinite(parameter_value) or parameter_value <= 0.0:
        raise ValueError(
            f'integ.{parameter_name} must be finite and positive, '
            f'got {parameter_value!r}'
        )

# seed points: circle on z=0 plane at r = 1.3*Rp
Npoints = 20
Rc = 1.3*Rp
theta = np.linspace(0,2*np.pi,Npoints,endpoint=False)
x0 = Rc*np.cos(theta)
y0 = Rc*np.sin(theta)
z0 = np.zeros_like(x0)
r0 = np.column_stack((x0,y0,z0))
Npoints = r0.shape[0]

# seed points: load from file
#file_input_seed_points = os.path.join(os.getenv('HOME'),'path/to/run/');
#r0 = np.loadtxt(file_input_seed_points,skiprows=0,usecols=(0, 1, 2),ndmin=2,)
#if r0.shape[1] != 3 or not np.all(np.isfinite(r0)):
#    raise ValueError(
#        f'Invalid seed-point file {file_input_seed_points}: expected three finite columns'
#    )
#Npoints = len(r0)

# square inner stopping radius
Rstop2 = Rstop * Rstop

# open file and check variable
vr = rhb.vlsvfile.VlsvReader(file_input_state)
if vr.check_variable(var_name) == False:
    print('ERROR: variable not found (' + var_name + ' / ' + file_input_state + ')')
    sys.exit(1)
var_dim = vr.read_variable_vectorsize(var_name)
if var_dim != 3:
    print('ERROR: variable should be a 3D vector (' + str(var_dim) + ' / ' + var_name + ')')
    sys.exit(1)

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
if not ( np.isclose(dx, dy, rtol=1e-12, atol=0.0) and np.isclose(dx, dz, rtol=1e-12, atol=0.0) ):
    print(f'ERROR: grid cells must be cubic (dx={dx/Rp}, dy={dy/Rp}, dz={dz/Rp})')
    sys.exit(1)
# cell center coordinates
x_centers = xmin + (np.arange(nx) + 0.5)*dx
y_centers = ymin + (np.arange(ny) + 0.5)*dx
z_centers = zmin + (np.arange(nz) + 0.5)*dx

# integrator parameters
coordinate_scale = max(
    abs(xmin), abs(xmax),abs(ymin),
    abs(ymax),abs(zmin),
    abs(zmax),
    dx,
)
target_error = integ.error_cells * dx
integ.max_step = integ.max_step_cells * dx
integ.atol = 0.5 * target_error
integ.rtol = 0.5 * target_error / coordinate_scale

# if prefer setting these manually, edit these
#integ.max_step = 7e3 # maximum step size [m]
#integ.atol = 1e-6 # absolute tolerance
#integ.rtol = 1e-3 # relative tolerance

# Read and reorder the vector field.  Vgrid uses (z, y, x, component) order;
# the final dimension lets one interpolation return all three components.
V = vr.read_variable(var_name)
cellids = vr.read_variable('CellID')

ncells = nx*ny*nz
if V.shape != (ncells, 3) or cellids.size != ncells:
    raise ValueError(
        f'Mesh/data size mismatch: expected {ncells} vector cells, '
        f'got V.shape={V.shape} and {cellids.size} CellIDs'
    )
cellids_sorted = np.argsort(cellids)
sorted_ids = cellids[cellids_sorted]
if np.any(np.diff(sorted_ids) != 1):
    raise ValueError('CellIDs must be unique and consecutive')
Vgrid = V[cellids_sorted, :].reshape(nz, ny, nx, 3)
del V, cellids, cellids_sorted, sorted_ids

# Vector-valued SciPy interpolator. This is also the reference implementation
# used to validate the compiled uniform-grid interpolator.
interpolator_scipy = sci_intp.RegularGridInterpolator(
    (z_centers, y_centers, x_centers),
    Vgrid,
    bounds_error=False,
    fill_value=np.nan,
    method=integ.intpol_method,
)


def _trilinear_uniform_python(field, x, y, z, x0, y0, z0, spacing):
    """Trilinear vector interpolation on a cell-centred uniform grid."""
    nz_local, ny_local, nx_local, ncomp = field.shape
    gx = (x - x0) / spacing
    gy = (y - y0) / spacing
    gz = (z - z0) / spacing

    if (gx < 0.0 or gx > nx_local - 1 or
            gy < 0.0 or gy > ny_local - 1 or
            gz < 0.0 or gz > nz_local - 1):
        return np.full(ncomp, np.nan)

    ix = int(np.floor(gx))
    iy = int(np.floor(gy))
    iz = int(np.floor(gz))
    if ix == nx_local - 1:
        ix = nx_local - 2
        tx = 1.0
    else:
        tx = gx - ix
    if iy == ny_local - 1:
        iy = ny_local - 2
        ty = 1.0
    else:
        ty = gy - iy
    if iz == nz_local - 1:
        iz = nz_local - 2
        tz = 1.0
    else:
        tz = gz - iz

    result = np.empty(ncomp, dtype=np.float64)
    for component in range(ncomp):
        c000 = field[iz,     iy,     ix,     component]
        c001 = field[iz,     iy,     ix + 1, component]
        c010 = field[iz,     iy + 1, ix,     component]
        c011 = field[iz,     iy + 1, ix + 1, component]
        c100 = field[iz + 1, iy,     ix,     component]
        c101 = field[iz + 1, iy,     ix + 1, component]
        c110 = field[iz + 1, iy + 1, ix,     component]
        c111 = field[iz + 1, iy + 1, ix + 1, component]

        c00 = c000 + tx * (c001 - c000)
        c01 = c010 + tx * (c011 - c010)
        c10 = c100 + tx * (c101 - c100)
        c11 = c110 + tx * (c111 - c110)
        c0 = c00 + ty * (c01 - c00)
        c1 = c10 + ty * (c11 - c10)
        result[component] = c0 + tz * (c1 - c0)
    return result


if njit is not None:
    trilinear_uniform_numba = njit(cache=True)(_trilinear_uniform_python)
else:
    trilinear_uniform_numba = None


def _unit_direction_uniform_python(field, x, y, z, x0, y0, z0, spacing):
    """Interpolate and normalize inside one compiled call for the RK45 RHS."""
    vector = trilinear_uniform_numba(
        field, x, y, z, x0, y0, z0, spacing
    )
    magnitude_squared = (
        vector[0] * vector[0]
        + vector[1] * vector[1]
        + vector[2] * vector[2]
    )
    if not np.isfinite(magnitude_squared):
        raise ValueError('Field interpolation is undefined')
    if magnitude_squared == 0.0:
        raise ValueError('Field direction is undefined at zero field')

    inverse_magnitude = 1.0 / np.sqrt(magnitude_squared)
    vector[0] *= inverse_magnitude
    vector[1] *= inverse_magnitude
    vector[2] *= inverse_magnitude
    return vector


if njit is not None:
    unit_direction_uniform_numba = njit(cache=True)(
        _unit_direction_uniform_python
    )
else:
    unit_direction_uniform_numba = None


def interpolate_scipy(crds):
    """Return the unnormalized SciPy-interpolated field vector."""
    x, y, z = crds
    return np.asarray(interpolator_scipy((z, y, x)), dtype=np.float64)


def interpolate_numba(crds):
    """Return the unnormalized compiled-interpolator field vector."""
    if trilinear_uniform_numba is None:
        raise RuntimeError(
            'The numba interpolation backend was selected, but numba is not installed. '
            "Install numba or set integ.intpol_backend = 'scipy'."
        )
    x, y, z = crds
    return trilinear_uniform_numba(
        Vgrid, x, y, z,
        x_centers[0], y_centers[0], z_centers[0], dx,
    )


def direction_numba(crds):
    """Return a unit field vector with interpolation and normalization compiled."""
    x, y, z = crds
    return unit_direction_uniform_numba(
        Vgrid, x, y, z,
        x_centers[0], y_centers[0], z_centers[0], dx,
    )


def _unit_direction(vector, crds):
    """Normalize the SciPy reference vector and check for invalid fields."""
    magnitude = np.sqrt(np.dot(vector, vector))
    if not np.isfinite(magnitude):
        raise ValueError(f'Field interpolation is undefined at {crds/Rp}')
    if magnitude == 0:
        raise ValueError(f'Field direction is undefined at zero field: {crds/Rp}')
    return vector / magnitude


def field_direction_scipy(t, crds):
    return _unit_direction(interpolate_scipy(crds), crds)


def field_direction_numba(t, crds):
    # solve_ivp requires the (t, state) signature; arc length t is unused.
    return direction_numba(crds)


if integ.intpol_backend == 'numba':
    if trilinear_uniform_numba is None:
        raise RuntimeError('Numba is not available')
    # Compile before tracing so the one-time JIT cost is excluded from tracing.
    # A synthetic nonzero field avoids making startup depend on the value in an
    # arbitrary simulation cell. Matching Vgrid's dtype also compiles the same
    # specialization that tracing will use.
    warmup_field = np.ones((2, 2, 2, 3), dtype=Vgrid.dtype)
    unit_direction_uniform_numba(
        warmup_field, 0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 1.0
    )
    del warmup_field
    field_direction = field_direction_numba
else:
    field_direction = field_direction_scipy

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


def trace_one_seed(seed, direction_function):
    """Trace one seed in both directions and return its line, class, and stats.

    The returned line is ordered from the backward endpoint, through the seed,
    to the forward endpoint.  Connectivity counts which directions hit Rstop.
    """
    start = perf_counter()
    stats = {
        'status': 'ok',
        'message': '',
        'forward_steps': 0,
        'backward_steps': 0,
        'forward_nfev': 0,
        'backward_nfev': 0,
        'elapsed': 0.0,
    }

    if np.linalg.norm(seed) <= Rstop:
        stats['status'] = 'skipped'
        stats['message'] = f'seed is at or inside Rstop ({Rstop/Rp:.3f} Rp)'
        stats['elapsed'] = perf_counter() - start
        return None, 3, stats
    if event_outer_boundaries(0, seed) <= 0:
        stats['status'] = 'skipped'
        stats['message'] = f'seed is outside the tracing box ({seed/Rp})'
        stats['elapsed'] = perf_counter() - start
        return None, 3, stats

    seed_as_float = seed.astype(float)
    sol_forward = None
    sol_backward = None
    try:
        sol_forward = sci_int.solve_ivp(
            direction_function,
            [0, integ.max_length],
            seed_as_float,
            method=integ.method,
            rtol=integ.rtol,
            atol=integ.atol,
            max_step=integ.max_step,
            events=[event_inner_boundary, event_outer_boundaries],
            vectorized=False,
        )
        sol_backward = sci_int.solve_ivp(
            lambda t, Y: -direction_function(t, Y),
            [0, integ.max_length],
            seed_as_float,
            method=integ.method,
            rtol=integ.rtol,
            atol=integ.atol,
            max_step=integ.max_step,
            events=[event_inner_boundary, event_outer_boundaries],
            vectorized=False,
        )
    except ValueError as exc:
        stats['status'] = 'failed'
        stats['message'] = str(exc)
        if sol_forward is not None:
            stats['forward_steps'] = len(sol_forward.t)
            stats['forward_nfev'] = sol_forward.nfev
        stats['elapsed'] = perf_counter() - start
        return None, 3, stats

    stats['forward_steps'] = len(sol_forward.t)
    stats['backward_steps'] = len(sol_backward.t)
    stats['forward_nfev'] = sol_forward.nfev
    stats['backward_nfev'] = sol_backward.nfev

    if not sol_forward.success or not sol_backward.success:
        stats['status'] = 'failed'
        stats['message'] = (
            f'forward={sol_forward.message}; backward={sol_backward.message}'
        )
        stats['elapsed'] = perf_counter() - start
        return None, 3, stats

    # Event index zero is the inner spherical boundary. Event index one is the
    # outer box boundary. If neither event occurs, the branch merely exhausted
    # max_length, so its eventual connectivity is unknown.
    hit_fwd = sol_forward.t_events[0].size > 0
    hit_bwd = sol_backward.t_events[0].size > 0
    hit_outer_fwd = sol_forward.t_events[1].size > 0
    hit_outer_bwd = sol_backward.t_events[1].size > 0
    limit_fwd = not hit_fwd and not hit_outer_fwd
    limit_bwd = not hit_bwd and not hit_outer_bwd
    n_hits = int(hit_fwd) + int(hit_bwd)
    if limit_fwd or limit_bwd:
        # A length-limited branch is inconclusive: it might reach either
        # boundary if tracing were continued. Do not label it as disconnected.
        connection = 3
        limited_directions = []
        if limit_fwd:
            limited_directions.append('forward')
        if limit_bwd:
            limited_directions.append('backward')
        branch_word = 'branches' if len(limited_directions) > 1 else 'branch'
        stats['status'] = 'length-limited'
        stats['message'] = (
            f'{" and ".join(limited_directions)} {branch_word} reached '
            f'max_length={integ.max_length/Rp:.3f} Rp without a boundary event'
        )
    elif n_hits == 2:
        connection = 2
    elif n_hits == 1:
        connection = 1
    else:
        connection = 0

    forward_trajectory = sol_forward.y.T
    backward_trajectory = sol_backward.y.T[::-1, :]
    field_line = np.vstack((backward_trajectory, forward_trajectory[1:, :]))
    stats['elapsed'] = perf_counter() - start
    return field_line, connection, stats


# Trace all seeds and render the result only when run as a script. Keeping this
# in main() lets the separate validation program import the numerical kernels.
def main():
    print('integrator parameters:')
    for name, value in vars(integ).items():
        print(name, '=', value)
    print('box_margin_cells =', box_margin_cells)
    print('tracing...')
    field_lines_all = []
    connectivity = []
    total_forward_steps = 0
    total_backward_steps = 0
    total_forward_nfev = 0
    total_backward_nfev = 0
    tracing_start = perf_counter()
    for ii, seed in enumerate(r0, start=1):
        print('field line ' + str(ii) + "/" + str(Npoints))
        field_line, conn, stats = trace_one_seed(seed, field_direction)
        field_lines_all.append(field_line)
        connectivity.append(conn)
        total_forward_steps += stats['forward_steps']
        total_backward_steps += stats['backward_steps']
        total_forward_nfev += stats['forward_nfev']
        total_backward_nfev += stats['backward_nfev']
        if stats['status'] != 'ok':
            print(
                f'  {stats["status"]} after {stats["elapsed"]:.3f} s: '
                f'{stats["message"]}'
            )
            continue
        print(
            f'  time={stats["elapsed"]:.3f} s, '
            f'nfev={stats["forward_nfev"] + stats["backward_nfev"]:,}, '
            f'accepted points={stats["forward_steps"] + stats["backward_steps"]:,}'
        )
    print(f'Tracing complete. Total forward steps: {total_forward_steps}, '
          f'total backward steps: {total_backward_steps}.')
    print(
        f'Total field evaluations: '
        f'{total_forward_nfev + total_backward_nfev:,} '
        f'(forward={total_forward_nfev:,}, backward={total_backward_nfev:,})'
    )
    print_timing('trace all field lines', tracing_start)

    # PLOTTING
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
        Line2D([0], [0], color='gray', lw=2,
               label='Invalid / incomplete / not traced'),
    ]
    ax.legend(handles=legend_elements, loc='best', title='Connectivity')
    plt.tight_layout()
    if save_output_figure == True:
        plt.savefig(file_output_figure)
        plt.clf()
        plt.close()
    else:
        plt.show(block=True)


if __name__ == '__main__':
    main()
