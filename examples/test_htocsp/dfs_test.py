"""
DFS sampler
"""

from __future__ import annotations
from time import time
import os
from pyxtal.optimize import DFS
from pyxtal.optimize.common import randomizer
from pyocse.lmp.pyxtal_calculator import LMP
from pyocse.parameters import ForceFieldParameters
exe = "/home/shattori/app/lammps_build/mylammps/build_fix_symmetry/lmp"

def lattice_cut_and_rotate(xtal, cut=2.0, ncut=3, verbose=False):
    if verbose:
        xtal.to_file('init.cif')
        print(xtal.get_separations([[1, 0, 0], [0, 1, 0], [0, 0, 1]]))
        print(xtal)
    for i in range(ncut): 
        xtal.cut_lattice(cut, True)
        if verbose:
            print(xtal.get_separations([[1, 0, 0], [0, 1, 0], [0, 0, 1]]))
            print(xtal)
            xtal.to_file(f'cut_{i}.cif')

    for site in xtal.mol_sites:
        site.optimize_orientation_by_energy(20, verbose=verbose)
    if verbose:
        xtal.to_file('opt.cif')
    return xtal

def myoptimizer(
    struc0,
    atom_info,
    workdir,
    tag = "job_0",
    opt_lat = True,
    max_time = 180,
    skip_ani = True,
):
    cwd = os.getcwd()
    t0 = time()
    os.chdir(workdir)

    v0 = struc0.lattice.volume
    struc_opt = lattice_cut_and_rotate(struc0, cut=2.0, ncut=3, verbose=False)
    vopt = struc_opt.lattice.volume
    struc = struc_opt.copy()
    chm_atom_info, lmp_atom_info = atom_info
    calc = LMP(struc, atom_info=lmp_atom_info, label=tag, exe=exe)
    calc.run(clean=False)

    results = {}
    results["xtal"] = calc.structure
    results["energy"] = calc.structure.energy
    results["time"] = time() - t0

    v_lmp = calc.structure.lattice.volume

    print("V0:", v0, "Vopt:", vopt, "LMP:", v_lmp)

    #comapare  with CHARMM
    #struc = struc0.copy()
    #t0 = time()
    #calc_chm = CHARMM(struc, steps=[1000,1000], atom_info=chm_atom_info, label=tag)
    #calc_chm.run(clean=False)
    #v_chm = calc_chm.structure.lattice.volume

    #print("V0:", v0, " LMP:", v_lmp, " CHARMM:", v_chm)
    #print("Eopt LMP:", results["energy"], " CHARMM:", calc_chm.structure.energy)
    #print("time LMP:", results["time"], " CHARMM:", time() - t0)

    os.chdir(cwd)

    return results




if __name__ == "__main__":
    import argparse
    import os

    from pyxtal.db import database

    calculator = LMP
    calculator.exe = exe

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-g",
        "--gen",
        dest="gen",
        type=int,
        default=5,
        help="Number of generation, optional",
    )
    parser.add_argument(
        "-p",
        "--pop",
        dest="pop",
        type=int,
        default=5,
        help="Population size, optional",
    )
    parser.add_argument("-n", "--ncpu", dest="ncpu", type=int,
                        default=1, help="cpu number, optional")
    parser.add_argument("--ffopt", action="store_true",
                        help="enable ff optimization")

    options = parser.parse_args()
    gen = options.gen
    pop = options.pop
    ncpu = options.ncpu
    ffopt = options.ffopt
    db_name, name = "test.db", "ACSALA" #please copy from pyxtal repo
    #db_name, name = "test.db", "DURNAH" #please copy from pyxtal repo
    #db_name, name = "test.db", "JAYDUI" #please copy from pyxtal repo

    wdir = name
    os.makedirs(wdir, exist_ok=True)
    os.makedirs(wdir + "/calc", exist_ok=True)

    db = database(db_name)
    row = db.get_row(name)
    xtal = db.get_pyxtal(name)
    smile, wt, spg = row.mol_smi, row.mol_weight, row.space_group.replace(
        " ", "")
    chm_info = None

    pmg0 = xtal.to_pymatgen()
    if xtal.has_special_site():
        xtal = xtal.to_subgroup()
    N_torsion = xtal.get_num_torsions()

    params = ForceFieldParameters([mol.smile for mol in xtal.molecules], style="gaff", chargemethod="am1bcc")

    print(name, xtal.group.number, xtal.has_special_site(), xtal.get_zprime())
    # GO run
    t0 = time()
    class myDFS (DFS):
        def _get_local_optimization_args(self):
            args = [
                randomizer,
                myoptimizer,
                self.smiles,
                self.block,
                self.num_block,
                (self.atom_info, self.direct_atominfo), #params instead of atom_info self.atom_info,
                self.workdir + "/" + "calc",
                self.sg,
                self.composition,
                self.lattice,
                self.torsions,
                self.molecules,
                self.sites,
                self.ref_pmg,
                self.matcher,
                self.ref_pxrd,
                self.use_hall,
                self.skip_ani,
                self.check_stable,
            ]
            return args
    #go = DFS(
    go = myDFS(
        smile,
        wdir,
        xtal.group.number,
        name.lower(),
        info=chm_info,
        ff_style="openff",  # 'gaff',
        ff_opt=ffopt,
        N_gen=gen,
        N_pop=pop,
        N_cpu=ncpu,
        cif="pyxtal.cif",
        skip_ani = True,
    )
    go.direct_atominfo = params

    suc_rate = go.run(pmg0)
    print(f"CSD {name:s} in Gen {go.generation:d}")

    if len(go.matches) > 0:
        best_rank = go.print_matches()
        mytag = f"True {best_rank:d}/{go.N_struc:d} Succ_rate: {suc_rate:7.4f}%"
    else:
        mytag = f"False 0/{go.N_struc:d}"

    eng = go.min_energy
    t1 = int((time() - t0)/60)
    strs = "Final {:8s} [{:2d}]{:10s} ".format(name, sum(xtal.numMols), spg)
    strs += "{:3d}m {:2d} {:6.1f}".format(t1, N_torsion, wt)
    strs += "{:12.3f} {:20s} {:s}".format(eng, mytag, smile)
    print(strs)
