import os
from ase import Atoms
#from ase.build import bulk
from ase.io import read
from ase.calculators.eam import EAM
from ase.calculators.emt import EMT
from ase.constraints import FixSymmetry, ExpCellFilter
from ase.optimize import LBFGS, FIRE
from ase.io import write
from ase.spacegroup.symmetrize import check_symmetry



eam_potential_file = 'Cu_u3.eam'
if not os.path.isfile(eam_potential_file):
    raise FileNotFoundError("eam_potential_file not found")

calc_eam = EAM(potential=eam_potential_file, elements=['Cu'])
calc_emt = EMT()

lat = 10.0
cu4 = Atoms(symbols=['Cu']*4, positions=[(0, 0, 0), (2.5, 0.0, 0.0), (2.5, 5.0, 0.0), (5.0, 5.0, 0.0)], pbc=True)
#cu4 = Atoms(symbols=['Cu']*4, positions=[(0, 0, 0), (0.76, 1.76, 0.55), (1.76, 0, 1.76), (0, 15.0, 15.0)], pbc=True)
cu4.cell = [[lat, 0, 0], [0, lat, 0], [0, 0, lat]]
# Generate bulk structure of Cu
#a0 = 3.615  # Actual lattice constant of Cu (Å)
#a0 = 3.75  # Actual lattice constant of Cu (Å)
#cubulk = bulk('Cu', 'fcc', a=a0, cubic=True)
#cubulk = cubulk.repeat([2, 2, 2])
cubulk = read('cu32.dat', format='lammps-data')
cubulk.symbols = ['Cu']*32
print(cubulk.cell.cellpar())
#cubulk.cell *= 1.1
# -------------------------------
# Main script
# -------------------------------


def custom_log_function(opt):
    step = opt.nsteps
    energy = opt.atoms.get_potential_energy()
    forces = opt.atoms.get_forces()
    print(f"Step: {step}, Energy: {energy:.4f} eV, Forces: {forces[0]}")

def save_xyz(atoms, step, fname):
    write(fname, atoms, append=True)
    #print(f"Step {step}: Structure saved to {filename}")

# Structure optimization (BFGS method)
def optimize_with_fix_symmetry(atoms, calc,
        fmax=1e-6, max_steps = 100,
        logfile = 'optimization.log',
        output_structure = 'optimized_structure.xyz'
    ):
    
    # Set calculator
    atoms.calc = calc
    
    # Apply FixSymmetry to maintain symmetry
    constraint = FixSymmetry(atoms, symprec=1e-4)
    atoms.set_constraint(constraint)
    ecf = ExpCellFilter(atoms)
    
    # remove logfile and output_structure if they exist
    if os.path.exists(logfile):
        os.remove(logfile)
    if os.path.exists(output_structure):
        os.remove(output_structure)
        
    # Structure optimization (BFGS method)
    opt = LBFGS(ecf, logfile=logfile)
    #opt = FIRE(atoms, a=0.1, logfile=logfile)
    opt.attach(lambda: save_xyz(atoms, opt.nsteps, output_structure), interval=1)
    opt.attach(lambda: custom_log_function(opt), interval=1)

    print("------\nOptimization start.")
    check_symmetry(atoms, symprec=1e-4, verbose=True)
    opt.run(fmax=fmax, steps=max_steps)
    
    # Save results
    #write(output_structure, atoms)
    
    # Display results
    print("------\nOptimization is complete.")
    check_symmetry(atoms, symprec=1e-4, verbose=True)
    raise
    print(f"Lattice constants and angles after optimization: {atoms.cell.cellpar()}")
    print(f"Energy after optimization: {atoms.get_potential_energy():.4f} eV")
    #print(f"Optimized structure is saved in '{output_structure}'.")
    #print(f"Optimization log is saved in '{logfile}'.")

if __name__ == '__main__':
    #optimize_with_fix_symmetry(cu4, calc_eam)
    optimize_with_fix_symmetry(cubulk, calc_emt)

