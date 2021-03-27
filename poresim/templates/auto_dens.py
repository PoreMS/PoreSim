import os
import porems as pms
import poresim as ps
import poreana as pa
import matplotlib.pyplot as plt


if __name__ == "__main__":
    # Add Todos
    print("Finish fill scripts ...")
    print("Finish ana.sh file ...")
    print("Add following script to running shell")
    # cd ../ana
    # sh ana.sh
    # python ana.py

    # Load molecule
    mol = pms.Molecule("molecule", "MOLSHORT", inp="MOLLINK")

    # Set analysis
    ana_list = {}
    ana_list["MOLSHORT"]  = {"traj": "traj.xtc",  "dens": True, "diff": False, "mol": mol, "atoms": []}

    # Run analysis
    for ana_name, ana_props in ana_list.items():
        sample = pa.Sample("../_gro/pore_system.obj", ana_props["traj"], ana_props["mol"], ana_props["atoms"], is_nojump=True)
        if ana_props["dens"]:
            sample.init_density("dens_"+ana_name+".obj")
        if ana_props["diff"]:
            sample.init_diffusion_bin("diff_"+ana_name+".obj", bin_num=35)
        sample.sample(is_parallel=True)

    # Calculate density
    dens = pa.density.calculate("dens_MOLSHORT.obj", target_dens=TARGETDENS)

    # Create plot
    pa.density.plot(dens)
    plt.gcf().suptitle(r"In: $\rho=$"+"%7.3f"%dens["in"][3]+r" kg m$^{-3}$, Out: $\rho=$"+"%7.3f"%dens["ex"][3]+r" kg m$^{-3}$")
    plt.savefig("density.pdf", format="pdf", dpi=1000)

    # Fill and rerun
    dens_ex = dens["ex"][3]
    num_diff = dens["diff"]
    if num_diff > 10:
        ps.utils.copy("../_fill/fillBackup.sh", "../_fill/fill.sh")
        ps.utils.replace("../_fill/fill.sh", "FILLDENS", str(int(num_diff)))

        os.system("cd ../_fill;sh fill.sh;cd ../min;SUBMITCOMMAND")
