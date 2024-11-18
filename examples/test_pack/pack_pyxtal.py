from pyxtal import pyxtal
from pyxtal.molecule import pyxtal_molecule
from pyocse.parameters import ForceFieldParameters
from pyxtal.interface.charmm import CHARMM
import subprocess
import lammps_logfile
import os

def lattice_cut_and_rotate(xtal, cut=2.0, ncut=10, verbose=False):
    xtal.to_file('init.cif')
    if verbose:
        print(xtal.get_separations([[1, 0, 0], [0, 1, 0], [0, 0, 1]]))
        print(xtal)
    for i in range(ncut):
        xtal.cut_lattice(cut, True)
        if verbose:
            print(xtal.get_separations([[1, 0, 0], [0, 1, 0], [0, 0, 1]]))
            print(xtal)

        for site in xtal.mol_sites:
            site.optimize_orientation_by_energy(verbose=verbose)
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
    asp = pyxtal_molecule('CC(=O)OC1=CC=CC=C1C(=O)O.smi', active_sites=[[11], [12],[20]])

    n_trial = 1
    for trial in range(n_trial):
        wdir = "Minimize_LMPSYM_random_ASP_" + str(trial)
        print("Start: ", wdir)
        os.makedirs(wdir, exist_ok=True)
        os.chdir(wdir)
        xtal0 = pyxtal(molecular=True)
        xtal0.from_random(3, 14, [asp], sites=[["4e"]])
        print(xtal0)

        xtal0 = lattice_cut_and_rotate(xtal0, cut=2.0, verbose=False)
        xtal = xtal0.copy()

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

        xtal = xtal0.copy()
        lmp_struc, _ = params.get_lmp_input_from_structure(xtal.to_ase(resort = False), xtal.numMols, set_template=False)
        lmp_struc.write_lammps(fin="pyxtal.in", fdat="pyxtal.dat")

        additional_lmpcmds_box = """
variable vmax equal 0.02
variable ptarget equal 25000
variable symprec equal 5e-4
min_style cg
min_modify dmax 0.02 line quadratic

fix br all symmetry ${symprec}
minimize 1e-5 1e-4 500 500 #1

unfix br
fix  br all box/relax/symmetry symprec ${symprec} x ${ptarget} y ${ptarget} z ${ptarget} xz 1 vmax ${vmax} fixedpoint 0 0 0 nreset 50
minimize 0 1e-4 500 500 #2

unfix br
fix  br all box/relax/symmetry symprec ${symprec} x 1 y ${ptarget} z 1 xz 1 vmax ${vmax} fixedpoint 0 0 0 nreset 50
minimize 0 1e-4 500 500 #3

unfix br
fix  br all box/relax/symmetry symprec ${symprec} x ${ptarget} y 1 z ${ptarget} xz ${ptarget} vmax ${vmax} fixedpoint 0 0 0 nreset 50
minimize 0 1e-4 500 500 #4

unfix br
fix br all box/relax/symmetry symprec  ${symprec} x 1 y 1 z 1 xz 1 vmax ${vmax} fixedpoint 0 0 0 nreset 50
minimize 0 1e-6 500 500 #5

unfix br
fix br all symmetry 1e-4
minimize 0 1e-6 200 200 #2
        """

        lmpintxt = open('pyxtal.in').read()
        lmpintxt = lmpintxt.replace("custom step ", "custom step spcpu ")
        lmpintxt = lmpintxt.replace("#compute ", "compute ")
        lmpintxt = lmpintxt.replace("#dump ", "dump ")
        lmpintxt = lmpintxt.replace("#dump_modify ", "dump_modify ")
        lmpintxt = lmpintxt.replace("kspace_modify ", "#kspace_modify ")
        lmpintxt += additional_lmpcmds_box
        open('pyxtal_mod.in', 'w').write(lmpintxt)
        cmd = '../../lmp -in pyxtal_mod.in -log lmp.log > /dev/null'
        _ = subprocess.getoutput(cmd)
        cputime = subprocess.getoutput("tail -n 1 lmp.log").split()[-1]

        step, eng, cell, sg = lammps_read("lmp.log")
        print("SG=", sg, "LMP:E=", eng, ",LAT=", cell, ", @@", step, ", @", cputime)

        os.chdir('..')

