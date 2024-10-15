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

    # Set box size (x,y,z)
    box = []
    if not box:
        print("Set box dimension to run ana.py")
        exit()

    # Set analysis
    ana_list = {}
    {% for mol in mols -%}
    ana_list["{{mol.name }}"] = {"traj": "traj_{{mol.name }}.xtc", "dens_box": True, "mc_trans": False, "mc": False, "mol": mol_dict["{{mol.name }}"], "atoms": []}
    {% endfor %}

    # Run analysis
    for ana_name, ana_props in ana_list.items():
        if ana_props["dens_box"]:
            sample = pa.Sample(box, ana_props["traj"], ana_props["mol"], ana_props["atoms"], [1 for x in ana_props["atoms"]])
            sample.init_density("dens_"+ana_name+"_box.obj")
            sample.sample(is_parallel=True)

        if ana_props["mc_trans"]:
            sample = pa.Sample(box, ana_props["traj"], ana_props["mol"], ana_props["atoms"], [1 for x in ana_props["atoms"]])
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
        num_diff = {}
        {% for mol in mols2 -%}
        {% if mol.area %}
        data = pa.utils.load("dens_{{mol.name }}_box.obj")
        {% for areas in mol.area %}
        # Define area of the molecule
        index_low = list(data["data"]["ex_width"]).index({{areas[0] }})
        index_up = list(data["data"]["ex_width"]).index({{areas[1] }})

        #Calculate denstiy and different number of molecules to target
        dens["{{mol.name }}"] = pa.density.bins("dens_{{mol.name }}_box.obj", target_dens={{mol.target_dens }}, area=[[0,1], [index_low,index_up]])
        num_diff["{{mol.name }}"] = dens["{{mol.name }}"]["diff"]

        # Write position file for the molecules in the specific area
        with open("../_gro/" + "position_{}_area{}.dat".format("{{mol.name }}",{{loop.index - 1}}), "w") as file_out:
            for i in range(int(num_diff["{{mol.name }}"]/2)):
                out_string = str({{mol.box[0] }}/2) + " "
                out_string += str({{mol.box[1] }}/2) + " "
                out_string += str(({{areas[1] }}+{{areas[0] }})/2) + "\n"
                file_out.write(out_string)
            file_out.close() 
        {% endfor %} 
        {% else %}
        #Calculate denstiy and different number of molecules to target
        dens["{{mol.name }}"] = pa.density.bins("dens_{{mol.name }}_box.obj", target_dens={{mol.target_dens }}, area=[[0,1], [0,150]])
        num_diff["{{mol.name }}"] = dens["{{mol.name }}"]["diff"]

        # Write position file for the molecules
        with open("../_gro/" + "position_{}.dat".format("{{mol.name }}"), "w") as file_out:
            for i in range(int(num_diff["{{mol.name }}"]/2)):
                out_string = str({{mol.box[0] }}/2) + " "
                out_string += str({{mol.box[1] }}/2) + " "
                out_string += str({{mol.box[2] }}/2) + "\n"
                file_out.write(out_string)
            file_out.close()
        {% endif %}
        {% endfor %} 

        # Adjust fill.sh and start simulation again
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
