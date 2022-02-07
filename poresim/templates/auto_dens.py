import os
import porems as pms
import poresim as ps
import poreana as pa


if __name__ == "__main__":
    # Todo
    print("Finish fill scripts ...")
    print("Finish ana.sh file ...")
    print("Add following script to running shell")
    # cd ../ana
    # sh ana.sh
    # python ana.py

    # Load molecule
    mol = pms.Molecule("molecule", "MOLSHORT", inp="MOLLINK")

    ## Silanol
    sioh = pms.Molecule("sioh", "SL")
    sioh.add("Si", [0, 0, 0], name="Si1")
    sioh.add("O", 0, r=1.3, name="O1")
    sioh.add("H", 1, r=1.0, name="H1")

    # Load pore
    pore = pa.utils.load("../_gro/pore.yml")

    # Set analysis
    ana_list = {}
    ana_list["MOLSHORT"] = {"traj": "traj.xtc", "dens": True, "dens_box": True, "diff": False, "mc_trans": False, "mc": False, "mol": mol, "atoms": []}
    # ana_list["sioh"] = {"traj": "traj_sioh.xtc", "dens": True, "dens_box": True, "diff": False, "mc_trans": False, "mc": False, "mol": sioh, "atoms": ["O1"]}

    # Run analysis
    for ana_name, ana_props in ana_list.items():
        if ana_props["dens_box"]:
            sample = pa.Sample(pore["dimensions"], ana_props["traj"], ana_props["mol"], ana_props["atoms"], [1 for x in ana_props["atoms"]])
            sample.init_density("dens_"+ana_name+"_box.obj", remove_pore_from_res=True)
            sample.sample(is_parallel=True)

        if ana_props["dens"] or ana_props["diff"]:
            sample = pa.Sample("../_gro/pore.yml", ana_props["traj"], ana_props["mol"], ana_props["atoms"], [1 for x in ana_props["atoms"]])
            if ana_props["dens"]:
                sample.init_density("dens_"+ana_name+".obj", remove_pore_from_res=True)
            if ana_props["diff"]:
                sample.init_diffusion_bin("diff_"+ana_name+".obj", bin_num=35)
            sample.sample(is_parallel=True)

        if ana_props["mc_trans"]:
            sample = pa.Sample(pore["dimensions"], ana_props["traj"], ana_props["mol"], ana_props["atoms"], [1 for x in ana_props["atoms"]])
            sample.init_diffusion_mc("diff_"+ana_name+"_trans.obj", [1, 2, 5, 10, 20, 30, 40, 50, 60, 70])
            sample.sample(is_parallel=True)

        if ana_props["mc"]:
            model = pa.CosineModel("diff_"+ana_name+"_trans.obj", 6, 10)
            pa.MC().run(model, "diff_"+ana_name+"_mc_cos.obj", nmc_eq=1000000, nmc=1000000)

    # Automation
    is_auto = True
    if is_auto:
        # Calculate density - area is given in bins
        dens = pa.density.bins("dens_MOLSHORT.obj", target_dens=TARGETDENS, area=[[10, 90], [10, 50]])

        # Fill and rerun
        num_diff = dens["diff"]
        if num_diff > 10:
            ps.utils.copy("../_fill/fillBackup.sh", "../_fill/fill.sh")
            ps.utils.replace("../_fill/fill.sh", "FILLDENS", str(int(num_diff)))

            os.system("cd ../_fill;sh fill.sh;cd ../min;SUBMITCOMMAND")
