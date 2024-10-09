import contextlib
import os
import shutil
import subprocess
import numpy as np
from pyocse.parameters import ForceFieldParameters
import lammps_logfile

def lammps_read(fname, last_thermo_pos=-3, sym_pos=-2):
    log = lammps_logfile.File(fname)
    step = log.get("Step")
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

class LMP:
    """
    A calculator to perform oragnic crystal structure optimization in CHARMM.

    Args:
        - struc: pyxtal.molecular_crystal.molecular_crystal object
        - label (str): label for this calculation
        - prefix (str): prefix of this calculation
        - exe (str): charmm executable
    """

    def __init__(
        self,
        struc,
        label="_",
        prefix="pyxtal",
        exe="lmp",
        timeout=300,
    ):
        self.errorE = 1e+5
        self.error = False
        # check charmm Executable
        #if shutil.which(exe) is None:
        #    raise BaseException(f"{exe} is not installed")
        #else:
        self.exe = exe
        self.timeout = timeout

        # Files IO
        self.prefix = prefix
        self.label = label
        self.inp = self.prefix + ".in"
        self.dat = self.prefix + ".dat"
        self.log = self.label + ".log"
        self.dump = "dump.lammpstrj"
        self.folder = "LMP_" + self.label

        # Structure Manipulation
        #struc.resort()
        self.structure = struc

    def write(self):
        xtal = self.structure
        params = ForceFieldParameters([mol.smile for mol in xtal.molecules], style="gaff", chargemethod="am1bcc")
        lmp_struc, _ = params.get_lmp_input_from_structure(xtal.to_ase(resort=False), xtal.numMols)
        lmp_struc.write_lammps(self.inp, self.dat)

        additional_lmpcmds_box = """
min_style cg
fix 1 all symmetry 5e-5 false false true true
minimize 1e-5 1e-5 20 20

unfix 1
fix 2 all box/relax/symmetry symprec 5e-5 aniso 0.0001 vmax 0.002
minimize 1e-5 1e-5 500 500
unfix 2
fix 2 all box/relax/symmetry symprec 5e-5 tri 0.0001 vmax 0.0001
minimize 1e-6 1e-6 500 500
        """
        lmpintxt = open(self.inp).read()
        lmpintxt = lmpintxt.replace("lmp.dat", self.dat)
        lmpintxt = lmpintxt.replace("custom step ", "custom step spcpu ")
        lmpintxt = lmpintxt.replace("#compute ", "compute ")
        lmpintxt = lmpintxt.replace("#dump ", "dump ")
        lmpintxt = lmpintxt.replace("#dump_modify ", "dump_modify ")
        lmpintxt += additional_lmpcmds_box
        open(self.inp, 'w').write(lmpintxt)
    
    def read(self):
        from ase.io import read
        step, eng, cell, sg = lammps_read(self.log)
        self.energy = eng
        ase_struc = read(self.dump, format='lammps-dump', index=-1)
        positions = ase_struc.get_positions()

        count = 0
        for _i, site in enumerate(self.structure.mol_sites):
            coords = positions[count: count + len(site.molecule.mol)]
            site.update(coords, self.structure.lattice)
            count += len(site.molecule.mol)
        # print("after relaxation  : ", self.structure.lattice, "iter: ", self.structure.iter)
        self.structure.optimize_lattice()
        self.structure.update_wyckoffs()

    def run(self, clean=True):
        """
        Only run calc if it makes sense
        """
        if not self.error:
            os.makedirs(self.folder, exist_ok=True)
            cwd = os.getcwd()
            os.chdir(self.folder)

            self.write()  # ; print("write", time()-t0)
            res = self.execute()  # ; print("exe", time()-t0)
            if res is not None:
                self.read()  # ; print("read", self.structure.energy)
            else:
                self.structure.energy = self.errorE
                self.error = True
            if clean:
                self.clean()

            os.chdir(cwd)

    def execute(self):
        cmd = '{self.exe} -in {self.inp} -log {self.log} > /dev/null'
        # os.system(cmd)
        with open(os.devnull, 'w') as devnull:
            try:
                # Run the external command with a timeout
                result = subprocess.run(
                    cmd, shell=True, timeout=self.timeout, check=True, stderr=devnull)
                return result.returncode  # Or handle the result as needed
            except subprocess.CalledProcessError as e:
                print(f"Command '{cmd}' failed with return code {e.returncode}.")
                return None
            except subprocess.TimeoutExpired:
                print(f"External command {cmd} timed out.")
                return None

    def clean(self):
        os.remove(self.inp) if os.path.exists(self.inp) else None
        os.remove(self.dat) if os.path.exists(self.dat) else None
        os.remove(self.log) if os.path.exists(self.log) else None
        os.remove(self.dump) if os.path.exists(self.dump) else None
