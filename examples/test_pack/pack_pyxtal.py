from pyxtal import pyxtal
from pyxtal.molecule import pyxtal_molecule
from pyocse.parameters import ForceFieldParameters
from pyxtal.interface.charmm import CHARMM
import subprocess
import lammps_logfile
import os

def lattice_cut_and_rotate(xtal, cut=2.0, ncut=3, verbose=False):
    xtal.to_file('init.cif')
    if verbose:
        print(xtal.get_separations([[1, 0, 0], [0, 1, 0], [0, 0, 1]]))
        print(xtal)
    for i in range(ncut): 
        xtal.cut_lattice(cut, True)
        if verbose:
            print(xtal.get_separations([[1, 0, 0], [0, 1, 0], [0, 0, 1]]))
            print(xtal)
    xtal.to_file('cut.cif')

    for site in xtal.mol_sites:
        site.optimize_orientation_by_energy(20, verbose=verbose)
    xtal.to_file('opt.cif')
    return xtal

def lammps_read(fname, sym_pos=-1):

    log = lammps_logfile.File(fname)
    step = log.get("Step")
    if "[Sym" not in step[sym_pos]:
        for i in range(1, 10):
            if "[Sym" in step[sym_pos-i]:
                sym_pos = sym_pos - i
                break
    last_thermo_pos = sym_pos - 1

    spcpu = log.get("S/CPU")
    pe = log.get("PotEng")
    a = log.get("Cella")
    b = log.get("Cellb")
    c = log.get("Cellc")
    alp = log.get("CellAlpha")
    bet = log.get("CellBeta")
    gam = log.get("CellGamma")

    ret = [
        float(step[last_thermo_pos]) / float(spcpu[last_thermo_pos]),
        pe[last_thermo_pos],
        (" ").join([
            str(a[last_thermo_pos]),
            str(b[last_thermo_pos]),
            str(c[last_thermo_pos]),
            str(alp[last_thermo_pos]),
            str(bet[last_thermo_pos]),
            str(gam[last_thermo_pos]),
        ]),
        spcpu[sym_pos]
    ]
    return ret

if __name__ == "__main__":
    style = 'openff'
    chargemethod = 'am1bcc'
    asp = pyxtal_molecule('CC(=O)OC1=CC=CC=C1C(=O)O.smi', active_sites=[(11), (12)])

    n_trial = 10
    for trial in range(n_trial):
        wdir = "Minimize_LMPSYM_random_ASP_" + str(trial)
        print("Start: ", wdir)
        os.makedirs(wdir, exist_ok=True)
        os.chdir(wdir)
        xtal = pyxtal(molecular=True)
        xtal.from_random(3, 14, [asp], sites=[["4e"]])
        xtal = lattice_cut_and_rotate(xtal, cut=2.0, verbose=False)
        #xtal = lattice_cut_and_rotate(xtal, cut=2.0, verbose=True)

        if xtal.has_special_site():
            xtal = xtal.to_subgroup()
        params = ForceFieldParameters([mol.smile for mol in xtal.molecules], style=style, chargemethod=chargemethod)
        params0 = params.params_init.copy()
        chm_struc = params.get_ase_charmm(params0)
        chm_info = chm_struc.get_charmm_info()
        with open("pyxtal_gen.prm", "w") as prm:
            prm.write(chm_info["prm"])
        with open("pyxtal_gen.rtf", "w") as rtf:
            rtf.write(chm_info["rtf"])

        chm = CHARMM(xtal, steps=[2000, 2000], prefix="pyxtal_gen", atom_info=chm_info)
        chm.run(clean=False)
        print("SG=",chm.structure.group.number, ",CHM:E=", chm.structure.energy, ",LAT=", chm.structure.lattice, ", @", chm.cputime)

        lmp_struc, _ = params.get_lmp_input_from_structure(xtal.to_ase(resort = False), xtal.numMols, set_template=False)
        lmp_struc.write_lammps(fin="pyxtal.in", fdat="pyxtal.dat")

        additional_lmpcmds = """
fix 1 all symmetry 5e-5
min_style       cg
minimize 1e-5 0 20 20

fix 2 all box/relax aniso 0.0001 vmax 0.002
minimize 0 1e-5 500 500
unfix 2
fix 2 all box/relax tri 0.0001 vmax 0.0001
minimize 1e-6 1e-6 500 500
        """
        additional_lmpcmds_box = """
min_style cg
fix 1 all symmetry 5e-5
minimize 1e-5 1e-5 20 20

unfix 2
fix 2 all box/relax/symmetry symprec 5e-5 aniso 0.0001 vmax 0.002
minimize 1e-5 1e-5 500 500
unfix 2
fix 2 all box/relax/symmetry symprec 5e-5 tri 0.0001 vmax 0.0001
minimize 1e-6 1e-6 500 500
        """
        lmpintxt = open('pyxtal.in').read()
        lmpintxt = lmpintxt.replace("custom step ", "custom step spcpu ")
        lmpintxt = lmpintxt.replace("#compute ", "compute ")
        lmpintxt = lmpintxt.replace("#dump ", "dump ")
        lmpintxt = lmpintxt.replace("#dump_modify ", "dump_modify ")
        lmpintxt += additional_lmpcmds_box
        open('pyxtal_mod.in', 'a').write(lmpintxt)
        cmd = '../../lmp -in pyxtal_mod.in -log lmp.log > /dev/null'
        _ = subprocess.getoutput(cmd)
        cputime = subprocess.getoutput("tail -n 1 lmp.log").split()[-1]

        step, eng, cell, sg = lammps_read("lmp.log")
        print("SG=", sg, "LMP:E=", eng, ",LAT=", cell, ", @@", step, ", @", cputime)

        os.chdir('..')
