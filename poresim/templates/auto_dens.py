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
    mol_dict = {}
    {% for mol in mols -%}
    mol_dict["{{mol.name }}"] = pms.Molecule("molecule", "{{mol.name }}", inp="{{mol.link }}")
    {% endfor %}


    ## Silanol
    sioh = pms.Molecule("sioh", "SL")
    sioh.add("Si", [0, 0, 0], name="Si1")
    sioh.add("O", 0, r=1.3, name="O1")
    sioh.add("H", 1, r=1.0, name="H1")

    # Load pore
    pore = pa.utils.load("../_gro/pore.yml")

    # Set analysis
    ana_list = {}
    {% for mol in mols -%}
    ana_list["{{mol.name }}"] = {"traj": "traj_{{mol.name }}.xtc", "dens": True, "dens_box": True, "diff": False, "mc_trans": False, "mc": False, "mol": mol_dict["{{mol.name }}"], "atoms": []}
    {% endfor %}

    # ana_list["sioh"] = {"traj": "traj_sioh.xtc", "dens": True, "dens_box": True, "diff": False, "mc_trans": False, "mc": False, "mol": sioh, "atoms": ["O1"]}

    # Run analysis
    for ana_name, ana_props in ana_list.items():
        if ana_props["dens_box"]:
            sample = pa.Sample(pore["system"]["dimensions"], ana_props["traj"], ana_props["mol"], ana_props["atoms"], [1 for x in ana_props["atoms"]])
            sample.init_density("dens_"+ana_name+"_box.obj", remove_pore_from_res=True)
            sample.sample(is_parallel=True)

        if ana_props["dens"] or ana_props["diff"]:
            sample = pa.Sample("../_gro/pore.yml", ana_props["traj"], ana_props["mol"], ana_props["atoms"], [1 for x in ana_props["atoms"]])
            if ana_props["dens"]:
                sample.init_density("dens_"+ana_name+".obj", remove_pore_from_res=False)
            if ana_props["diff"]:
                sample.init_diffusion_bin("diff_"+ana_name+".obj", bin_num=35)
            sample.sample(is_parallel=True)

        if ana_props["mc_trans"]:
            sample = pa.Sample("../_gro/pore.yml", ana_props["traj"], ana_props["mol"], ana_props["atoms"], [1 for x in ana_props["atoms"]])
            sample.init_diffusion_mc("diff_"+ana_name+"_trans.obj", [1, 2, 5, 10, 20, 30, 40, 50, 60, 70])
            sample.sample(is_parallel=True)

        if ana_props["mc"]:
            model = pa.CosineModel("diff_"+ana_name+"_trans.obj", 6, 10)
            pa.MC().run(model, "diff_"+ana_name+"_mc_cos.obj", nmc_eq=1000000, nmc=1000000)

    {% if fill %}
    # Automation
    is_auto = True
    if is_auto:
        # Calculate density - area is given in bins
        dens = {}
        {% for mol in mols2 -%}
        dens["{{mol.name }}"] = pa.density.bins("dens_{{mol.name }}.obj", target_dens={{mol.target_dens }}, area=[[0, 100], [0, 100]])
        {% endfor %}
        # Fill and rerun
        num_diff = {}
        {% for mol in mols2 -%}
        num_diff["{{mol.name }}"] = dens["{{mol.name }}"]["diff"]       
        with open("../_gro/" + "position_{}.dat".format("{{mol.name }}"), "w") as file_out:
            for i in range(int(num_diff["{{mol.name }}"]/2)):
                out_string = str(pore["system"]["dimensions"][0]/2) + " "
                out_string += str(pore["system"]["dimensions"][1]/2) + " "
                out_string += str(pore["system"]["reservoir"]/2) + "\n"
                file_out.write(out_string)
            for i in range(int(num_diff["{{mol.name }}"]/2)):
                out_string = str(pore["system"]["dimensions"][0]/2) + " "
                out_string += str(pore["system"]["dimensions"][1]/2) + " "
                out_string += str(pore["system"]["dimensions"][2] -pore["system"]["reservoir"]/2) + "\n"
                file_out.write(out_string)
            file_out.close()
        {% endfor %}
        if (all(i<10 for i in num_diff.values()))==False:
            ps.utils.copy("../_fill/fillBackup.sh", "../_fill/fill.sh")
            {% for mol in mols2 -%}
            if num_diff["{{mol.name }}"] > 10:
                ps.utils.replace("../_fill/fill.sh", "FILLDENS_{{mol.name }}", str(int(num_diff["{{mol.name }}"])))
            else:
                ps.utils.replace("../_fill/fill.sh", "FILLDENS_{{mol.name }}", str(int(0)))
            {% endfor %}
            os.system("cd ../_fill;sh fill.sh;cd ../min;{{submit }}")

    {% else %}
    # Nothing to fill up
    {% endif %}   
