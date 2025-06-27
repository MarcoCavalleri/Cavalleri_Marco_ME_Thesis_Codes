import pde
import numpy as np
from pde import PDE, FieldCollection, PlotTracker, ScalarField, UnitGrid, MemoryStorage, CartesianGrid # <-- Zwicker py-pde library
from numpy import exp, log, pi, savetxt, dot, gradient
import os



# === PARAMETERS ===

# Same as in the MatLab code

D = 1
D_1, D_2, D_3, D_s = D, D, D, D

mu0_1, mu0_2, mu0_3, mu0_s = 0, -2, 0, 0

s = 0.001
s_1, s_2, s_3, s_4, s_5 = s, s, s, s, s

# L and K ate not defined as matrices but it will be used their entries directly in the equations
chi = 1
k1 = 0.5
k2 = 0.1
mu_Y1, mu_Y2 = 1, -1

# Initial concentrations (homogeneous steady state)
c0_1 = 0.121606
c0_2 = 2.715039
c0_3 = 0.102156
c0_s = 1.000000

# Simulation parameters
# L is the size of the grid, n_L is the number of points in each direction
L = 30
n_L = 3*L

# Perturbation
ic_type = "random" # <-- Random or Previous
perturbation = 0.01 # <-- Square root of the variance

# Simulation time and storage interval
T = 1000
storage_interval = 100



# === PDE ===

# Chemical potentials
mu_1 = f"({mu0_1} + log(c_1) + {chi}*c_2-{k1}*laplace(c_1)-{k2}*laplace(c_2))"
mu_2 = f"({mu0_2} + log(c_2) + {chi}*(c_1+c_3+c_s) - {k1}*laplace(c_2) - {k2}*(laplace(c_1+c_3+c_s)))"
mu_3 = f"({mu0_3} + log(c_3) + {chi}*c_2-{k1}*laplace(c_3)-{k2}*laplace(c_2))"
mu_s = f"({mu0_s} + log(c_s) + {chi}*c_2-{k1}*laplace(c_s)-{k2}*laplace(c_2))"

# Reaction currents
j_1 = f"{s_1}*(exp({mu_Y1})-exp({mu_1}))"
j_2 = f"{s_2}*(exp({mu_1})-exp({mu_2}))"
j_3 = f"{s_3}*(exp({mu_2})-exp({mu_Y2}))"
j_4 = f"{s_4}*(exp({mu_2})-exp({mu_3}))"
j_5 = f"{s_5}*(exp({mu_3})-exp({mu_1}))"

# Chemical reaction terms
chem_term_1 = f"({j_1}-{j_2}+{j_5})"
chem_term_2 = f"({j_2}-{j_3}-{j_4})"
chem_term_3 = f"({j_4}-{j_5})"

# Diffusion terms
diff_term_1 = f"{D_1}*(dot(gradient(c_1),gradient({mu_1}))+c_1*laplace({mu_1}))"
diff_term_2 = f"{D_2}*(dot(gradient(c_2),gradient({mu_2}))+c_2*laplace({mu_2}))"
diff_term_3 = f"{D_3}*(dot(gradient(c_3),gradient({mu_3}))+c_3*laplace({mu_3}))"
diff_term_s = f"{D_s}*(dot(gradient(c_s),gradient({mu_s}))+c_s*laplace({mu_s}))"

# Define the PDE
# "c_i" represents the temporal derivative 
eq = PDE(
    {
        "c_1": f"{diff_term_1}+{chem_term_1}",
        "c_2": f"{diff_term_2}+{chem_term_2}",
        "c_3": f"{diff_term_3}+{chem_term_3}",
        "c_s": f"{diff_term_s}"
    }
)



# === INITIAL CONDITIONS ===

# Grid & PBCs
grid = CartesianGrid([[0,L], [0,L]], [n_L,n_L], periodic=[True, True])

physics = f"model-0-fig1_{ic_type}_T={T}"
name_of_calculation = f"{physics}_L={L}_Ngrid={n_L}_icPert={perturbation}"
print(name_of_calculation)

if ic_type == "previous":
    # Load initial conditions from files
    input_dir = r'results_fig1_ss2' # <-- Directory
    filenames = {
        'c_1': os.path.join(input_dir, 'c_1_fig1.txt'),
        'c_2': os.path.join(input_dir, 'c_2_fig1.txt'),
        'c_3': os.path.join(input_dir, 'c_3_fig1.txt'),
        'c_s': os.path.join(input_dir, 'c_s_fig1.txt'),
    }

    # Helper function to load concentration from file
    def load_concentration(filename, n_x, n_y):
        data = np.loadtxt(filename, comments='#')
        field_flat = data[:,2]
        field = field_flat.reshape((n_x, n_y))
        return field

    c_1_init_data = load_concentration(filenames['c_1'], n_L, n_L)
    c_2_init_data = load_concentration(filenames['c_2'], n_L, n_L)
    c_3_init_data = load_concentration(filenames['c_3'], n_L, n_L)
    c_s_init_data = load_concentration(filenames['c_s'], n_L, n_L)

    c_1 = ScalarField(grid, data=c_1_init_data, label="$c_1$")
    c_2 = ScalarField(grid, data=c_2_init_data, label="$c_2$")
    c_3 = ScalarField(grid, data=c_3_init_data, label="$c_3$")
    c_s = ScalarField(grid, data=c_s_init_data, label="$c_s$")
elif ic_type == "random":
    # Radom initial conditions
    c_1 = c0_1 * (1 + perturbation * ScalarField.random_normal(grid, label="$c_1$"))
    c_2 = c0_2 * (1 + perturbation * ScalarField.random_normal(grid, label="$c_2$"))
    c_3 = c0_3 * (1 + perturbation * ScalarField.random_normal(grid, label="$c_3$"))
    c_s = c0_s * (1 + perturbation * ScalarField.random_normal(grid, label="$c_s$"))

# Assemble
state = FieldCollection([c_1, c_2, c_3, c_s])



# === SIMULATION ===

storage = MemoryStorage()
printinterval = 1000
trackers = [
   "progress",  # show progress bar during simulation
   "steady_state",  # abort when steady state is reached
   storage.tracker(storage_interval),  # store data every simulation time unit
   PlotTracker(printinterval, plot_args = {"vmin": 0, "vmax": 1},  show = False ),  # show images during simulation
]

# Adaptive explicit Runge-Kutta scheme
solver = pde.ExplicitSolver(eq, scheme="runge-kutta", adaptive=True)
controller = pde.Controller(solver, t_range = T, tracker = trackers)
sol = controller.run(state, dt = 1e-7)
sol.label = "explicit rk solver"

print("Diagnostic information from second run:")
print(controller.diagnostics)

sol.plot()



# == SAVE RESULTS ===

output_dir = r"results_fig1_ss2" # <-- Possibly the same as input_dir
os.makedirs(output_dir, exist_ok=True)

# Final concentration fields
c_1_final = sol[0]
c_2_final = sol[1]
c_3_final = sol[2]
c_s_final = sol[3]

# Create the grid coordinates
x_coords = grid.axes_coords[0]
y_coords = grid.axes_coords[1]
X, Y = np.meshgrid(x_coords, y_coords, indexing='ij')

c_1_data = np.column_stack((X.ravel(), Y.ravel(), c_1_final.data.ravel()))
c_2_data = np.column_stack((X.ravel(), Y.ravel(), c_2_final.data.ravel()))
c_3_data = np.column_stack((X.ravel(), Y.ravel(), c_3_final.data.ravel()))
c_s_data = np.column_stack((X.ravel(), Y.ravel(), c_s_final.data.ravel()))

np.savetxt(os.path.join(output_dir, "c_1_fig1.txt"), c_1_data, header="x y c_1")
np.savetxt(os.path.join(output_dir, "c_2_fig1.txt"), c_2_data, header="x y c_2")
np.savetxt(os.path.join(output_dir, "c_3_fig1.txt"), c_3_data, header="x y c_3")
np.savetxt(os.path.join(output_dir, "c_s_fig1.txt"), c_s_data, header="x y c_s")
