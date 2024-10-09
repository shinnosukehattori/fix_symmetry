from pyocse.parameters import ForceFieldParameters
from pyxtal.interface.charmm import CHARMM
import subprocess
import os
import lammps_logfile

def lammps_read(fname, last_thermo_pos=-3, sym_pos=-2):
    log = lammps_logfile.File(fname)
    x = log.get("Step")
    y = log.get("PotEng")
    a = log.get("Cella")
    b = log.get("Cellb")
    c = log.get("Cellc")
    alp = log.get("CellAlpha")
    bet = log.get("CellBeta")
    gam = log.get("CellGamma")

    celpar_ary = [
        str(a[last_thermo_pos]),
        str(b[last_thermo_pos]),
        str(c[last_thermo_pos]),
        str(alp[last_thermo_pos]),
        str(bet[last_thermo_pos]),
        str(gam[last_thermo_pos]),
    ]

    return x[last_thermo_pos], y[last_thermo_pos], y[sym_pos], (" ").join(celpar_ary)

if __name__ == "__main__":
    from pyxtal.db import database

    # db = database('../HT-OCSP/benchmarks/Si.db')
    db = database("./test.db")
    style = 'gaff' #'openff'
    #chargemethod = 'gasteiger'
    chargemethod = 'am1bcc'
    #xtal = db.get_pyxtal("ACSALA")
    codes = db.get_all_codes()
    #codes = ["VOBYAN"]
    #codes = ["AXIDER"]
    #codes = ["ACSALA"]
    skippes = ["TIDFES", "XAFQAZ", "BOQQUT", "BOQQUT01", "UJIRIO01", "UJIRIO05"] #openff.toolkit.utils.exceptions.UndefinedStereochemistryError: Unable to make OFFMol from SMILES: RDMol has unspecified stereochemistry. Undefined chiral centers
    for code in codes:
        if code in skippes:
            continue
        #workdir
        wdir = "Minimize_LMPSYM_" + code
        os.makedirs(wdir, exist_ok=True)
        os.chdir(wdir)

        
        xtal = db.get_pyxtal(code)
        print(code, xtal.group.number, xtal.has_special_site(), xtal.get_zprime())

        # for charmm
        row = db.get_row(code)
        if xtal.has_special_site():
            xtal = xtal.to_subgroup()
        chm_info = row.data["charmm_info"]
        with open("pyxtal.prm", "w") as prm:
            prm.write(chm_info["prm"])
        with open("pyxtal.rtf", "w") as rtf:
                rtf.write(chm_info["rtf"])
        chm = CHARMM(xtal, steps=[2000], atom_info=chm_info)
        chm.run(clean=False)
        print(code, "SG=",chm.structure.group.number, ",CHM:E=", chm.structure.energy, ",LAT=", chm.structure.lattice, ", @", chm.cputime)

        # for lammps
        xtal = db.get_pyxtal(code)
        params = ForceFieldParameters([mol.smile for mol in xtal.molecules], style=style, chargemethod=chargemethod)
        lmp_struc, _ = params.get_lmp_input_from_structure(xtal.to_ase(resort=False), xtal.numMols)
        lmp_struc.write_lammps()

        additional_lmpcmds = """
fix 1 all symmetry 5e-5 false false true true true 
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
fix 1 all symmetry 5e-5 false false true true
minimize 1e-5 1e-5 20 20

unfix 2
fix 2 all box/relax/symmetry symprec 5e-5 aniso 0.0001 vmax 0.002
minimize 1e-5 1e-5 500 500
unfix 2
fix 2 all box/relax/symmetry symprec 5e-5 tri 0.0001 vmax 0.0001
minimize 1e-6 1e-6 500 500
        """
        lmpintxt = open('lmp.in').read()
        lmpintxt = lmpintxt.replace("#compute ", "compute ")
        lmpintxt = lmpintxt.replace("#dump ", "dump ")
        lmpintxt = lmpintxt.replace("#dump_modify ", "dump_modify ")
        lmpintxt += additional_lmpcmds_box
        open('lmpbox.in', 'w').write(lmpintxt)
        cmd = '../../lmp -in lmpbox.in -log lmpbox.log > /dev/null'
        _ = subprocess.getoutput(cmd)
        cputime = subprocess.getoutput("tail -n 1 lmpbox.log").split()[-1]

        step, eng, sg, cell = lammps_read("lmpbox.log")
        print(code, "SG=", sg, "LMP:E=", eng, ",LAT=", cell, ", @@", step, ", @", cputime)


        os.chdir('..')
