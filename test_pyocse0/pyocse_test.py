from pyocse.parameters import ForceFieldParameters
import subprocess
import os

if __name__ == "__main__":
    from pyxtal.db import database

    # db = database('../HT-OCSP/benchmarks/Si.db')
    db = database("./test.db")
    style = 'gaff' #'openff'
    chargemethod = 'gasteiger'
    #chargemethod = 'am1bcc'
    #xtal = db.get_pyxtal("ACSALA")
    codes = db.get_all_codes()
    #codes = ["VOBYAN"]
    codes = ["AXIDER"]
    for code in codes:
        #workdir
        wdir = "Minimize_LMPSYM_" + code
        os.makedirs(wdir, exist_ok=True)
        os.chdir(wdir)
        
        xtal = db.get_pyxtal(code)
        if xtal.has_special_site():
            xtal = xtal.to_subgroup()

        print(code, xtal.group.number)
        smiles = [mol.smile for mol in xtal.molecules]
        assert smiles[0] is not None
        params = ForceFieldParameters(smiles, style=style, chargemethod=chargemethod)
        params0 = params.params_init.copy()

        lmp_struc, _ = params.get_lmp_input_from_structure(xtal.to_ase(resort=False), xtal.numMols)
        lmp_struc.write_lammps()

        additional_lmpcmds = """
fix   sym all symmetry 2e-5
minimize 1e-6 1e-6 10 10

#fix   2 all box/relax aniso 1e-5 vmax 0.0001
#minimize 1e-6 1e-6 500 500

#unfix 2
#fix      2 all box/relax tri 1e-5 vmax 0.00001
#minimize 0 1e-8 500 500
        """
        open('lmp.in', 'a').write(additional_lmpcmds)
        cmd = '../../lmp -in lmp.in -log lmp.log'
        _ = subprocess.getoutput(cmd)

        #grep  finale space from lmp.log
        cmd = "grep 'Space group number' lmp.log"
        #execute the command and get the output
        output = subprocess.getoutput(cmd)
        print(output)
        

        os.chdir('..')
