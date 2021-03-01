import os

import matplotlib.pyplot as plt

import porems as pms
import poreana as pa


# Add Todos
print("Finish fill scripts ...")
print("Finish ana.sh file ...")
print("Add following script to running shell")
return
# cd ../ana
# sh ana.sh
# python ana.py

# Load molecule
mol = pms.Molecule("molecule", "MOLSHORT", inp="MOLLINK")

# Sample density
pa.density.sample("../_gro/pore_system.obj", "traj.pdb", "traj.trr", "dens.obj", mol)

# Calculate density
dens = pa.density.calculate("dens.obj", target_dens=TARGETDENS)

# Create plot
pa.density.plot(dens)
plt.gcf().suptitle(r"In: $\rho=$"+"%7.3f"%dens["in"][3]+r" kg m$^{-3}$, Out: $\rho=$"+"%7.3f"%dens["out"][3]+r" kg m$^{-3}$")
plt.savefig("density.pdf", format="pdf", dpi=1000)

# Fill and rerun
dens_out = dens["out"][3]
num_diff = dens["diff"]
if num_diff > 10:
    psa.utils.copy("../_fill/fillBackup.sh", "../_fill/fill.sh")
    psa.utils.replace("../_fill/fill.sh", "FILLDENS", str(int(num_diff)))

    os.system("cd ../_fill;sh fill.sh;cd ../min;SUBMITCOMMAND")
