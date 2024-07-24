Create a two Phase system
=========================
.. code:: ipython3

    import poresim as ps

Set job inputs for your working station
---------------------------------------

Example Bash file:

::

   #!/bin/bash 
   #SBATCH --nodes=SIMULATIONNODES 
   #SBATCH --ntasks-per-node=SIMULATIONPROCS 
   #SBATCH --time=SIMULATIONTIME 
   #SBATCH --job-name=SIMULATIONLABEL 
   #SBATCH --error=SIMULATIONLABEL.e.%J 
   #SBATCH --output=SIMULATIONLABEL.o.%J 

   module load chem/gromacs/2019.6_cpu 

   COMMANDCHANGEDIR

   COMMANDGROMACS

Simulation inputs

.. code:: ipython3

    boxes = {}
    
    job_bw = {"min": {"file": "data/bash.sh", "nodes": 1, "np": 1, "wall": "00:30:00"},
              "nvt": {"file": "data/bash.sh", "nodes": 1, "np": 20,"wall": "20:00:00"},
              "npt": {"file": "data/bash.sh","nodes": 1, "np": 20, "runs": 15, "wall": "20:00:00"},
              "run": {"file": "data/bash.sh", "maxh": 8, "nodes": 1, "np": 20, "runs": 15, "wall": "8:00:00"}}
    
    param = {"min": {"file": "../Simulation_Inputs/2Phase/2phase_eq_min.mdp"},
                 "nvt": {"file": "../Simulation_Inputs/2Phase/2phase_eq_nvt.mdp", "param": {"NUMBEROFSTEPS":   10000000, "TEMPERATURE_VAL": 353.15}}, 
                 "npt": {"file": "../Simulation_Inputs/2Phase/2phase_eq_npt.mdp", "param": {"NUMBEROFSTEPS":   30000000, "TEMPERATURE_VAL": 353.15}}, 
                 "run": {"file": "../Simulation_Inputs/2Phase/2phase_prod_nvt.mdp", "param": {"NUMBEROFSTEPS":  200000000, "TEMPERATURE_VAL": 353.15}}} 
    
    
    cluster = {"address": "",                    
                      "directory": "xxxxx/il/",        
                      "queuing": {"add_np": True, "mpi": "mpirun -n", "shell": "bash.sh", "submit": "sbatch -p cpuonly"}}

Create a box object with the molecules contained in your system
---------------------------------------------------------------

.. code:: ipython3

    # Add simulation label
    pores = ps.Box("353_2phase")
    pores.set_label("353_2phase")

    # Import empty gro file with the dimensions of box
    pores.add_box("data/box_2phase.gro")
    
    # Add gro files of the molecules
    # IL Phase
    pores.add_mol("CAT", "../Molecules/cat.gro", inp=10, area=[[0,5],[15,20]], box = [8,8,20], kwargs_gmx={"-try":1000, "-scale":0.47})
    pores.add_mol("EDC", "../Molecules/edc.gro", inp=10, area = [[5,15]], box = [8,8,20], kwargs_gmx={"-try":1000, "-scale":0.47})
    pores.add_mol("IM", "../Molecules/im.gro", inp=1600, area = [[0,5],[15,20]], box = [8,8,20], kwargs_gmx={"-try":1000, "-scale":0.47})
    pores.add_mol("BF4", "../Molecules/bf4.gro", inp=1620, area = [[0,5],[15,20]], box = [8,8,20], kwargs_gmx={"-try":1000, "-scale":0.47})
    
    
    #Heptane Phase
    pores.add_mol("HEP", "../Molecules/1-heptane.gro", inp=2420, area = [[5,15]], box = [8,8,20],kwargs_gmx={"-try":1000, "-scale":0.47})
    
    # Add top files 
    pores.add_topol("../Topologies/forcefield_new_red.itp", "top")
    pores.add_topol("../Topologies/ffnonbonded_new_red.itp","top")
    pores.add_topol("../Topologies/ffbonded_new.itp","top")
    pores.add_topol("../Topologies/cat.itp", "top")
    pores.add_topol("../Topologies/bf4.itp", "top")
    pores.add_topol("../Topologies/im.itp", "top")
    pores.add_topol("../Topologies/1-heptane.itp", "top")
    pores.add_topol("../Topologies/educt_oplsaa.itp", "top")
    pores.add_topol("../Topologies/master_topol_2phase.top", "master")
    
    # Set job and parameters
    pores.set_job(job_bw)
    pores.set_param(param)


Create simulation objects and generate folder structure
--------------------------------------------------------

.. code:: ipython3

    sim = ps.Simulate("Simulations/2phase/353_2phase", pores)  # Single
    sim.set_cluster(cluster)
    sim.generate()
    
    


``Finished simulation folder - 353_2phase ...``


Image of the system after filling the box
-----------------------------------------

Inside view of the pore after running construct.sh and energy minimization. 
Colour code: catalyst, red; ionic anion, orange; ionic kation, blue; rectant, grey; heptane, pink.

.. figure::  /pics/2phase_il.pdf
      :align: center
      :width: 50%
      :name: fig1
