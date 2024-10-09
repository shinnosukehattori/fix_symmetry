"""
WFS sampler
"""

from __future__ import annotations
from time import time
from pyxtal.optimize import WFS


if __name__ == "__main__":
    import argparse
    import os

    from pyxtal.db import database

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-g",
        "--gen",
        dest="gen",
        type=int,
        default=1,
        help="Number of generation, optional",
    )
    parser.add_argument(
        "-p",
        "--pop",
        dest="pop",
        type=int,
        default=1,
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
    wdir = name
    os.makedirs(wdir, exist_ok=True)
    os.makedirs(wdir + "/calc", exist_ok=True)

    db = database(db_name)
    row = db.get_row(name)
    xtal = db.get_pyxtal(name)
    smile, wt, spg = row.mol_smi, row.mol_weight, row.space_group.replace(
        " ", "")
    chm_info = None

    #if not ffopt:
    #    if "charmm_info" in row.data:
    #        # prepare charmm input
    #        chm_info = row.data["charmm_info"]
    #        with open(wdir + "/calc/pyxtal.prm", "w") as prm:
    #            prm.write(chm_info["prm"])
    #        with open(wdir + "/calc/pyxtal.rtf", "w") as rtf:
    #            rtf.write(chm_info["rtf"])
    #    else:
    #        # Make sure we generate the initial guess from ambertools
    #        if os.path.exists("parameters.xml"):
    #            os.remove("parameters.xml")
    # load reference xtal
    pmg0 = xtal.to_pymatgen()
    if xtal.has_special_site():
        xtal = xtal.to_subgroup()
    N_torsion = xtal.get_num_torsions()

    # GO run
    t0 = time()
    go = WFS(
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
    )

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
