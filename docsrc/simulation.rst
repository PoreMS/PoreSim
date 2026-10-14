:orphan:

.. raw:: html

  <div class="container-fluid">
    <div class="row">
      <div class="col-md-10">
        <div style="text-align: justify; text-justify: inter-word;">

Create a simulation box
=======================

.. code-block:: python

  import poresim as ps


Set job, parameter, and cluster dictionaries
--------------------------------------------

.. code-block:: python

  job = {"min": {"file": "data/horeka.sh", "nodes": 1, "np":  1, "wall": "00:30:00"},
         "nvt": {"file": "data/horeka.sh", "nodes": 2, "np": 75, "wall": "72:00:00"},
         "run": {"file": "data/horeka.sh", "maxh": 72, "nodes": 4, "np": 75, "runs": 5, "wall": "72:00:00"}}

  param = {"min": {"file": "data/mdp/pore_min.mdp"},
           "nvt": {"file": "data/mdp/pore_nvt.mdp", "param": {"NUMBEROFSTEPS":     20000, "TEMPERATURE_VAL": 295}},
           "run": {"file": "data/mdp/pore_run.mdp", "param": {"NUMBEROFSTEPS": 200000000, "TEMPERATURE_VAL": 295}}}

  cluster = {"address": "xx_xxxxxxxx@horeka.scc.kit.edu",
             "directory": "/link/to/simulation/folder/",
             "queuing": {"add_np": False, "mpi": "$DO_PARALLEL", "shell": "horeka.sh", "submit": "sbatch --partition cpuonly"}}


Create box objects
------------------

.. code-block:: python

  # Define molecules
  mols = {"methanol": "1OL", "ethanol": "2OL"}

  # Define target density
  dens = {}  # kg/m^3
  dens["1OL"] = 788.009
  dens["2OL"] = 792.485

   # Define molar mass
  mass = {}  # g/mol
  mass["1OL"] = 32.04
  mass["2OL"] = 46.068

  # Generate simulation boxes
  boxes = {}
  pore_link = "pore/"
  for mol, short in mols.items():
      sim_name = mol
      sim_label = short
      boxes[sim_label] = ps.Box(sim_name)
      boxes[sim_label].set_label(sim_label)
      boxes[sim_label].add_box(pore_link+"pore.gro")
      boxes[sim_label].add_pore(pore_link+"pore.yml")
      boxes[sim_label].add_struct("mol_list", pore_link+"pore.obj")
      boxes[sim_label].add_struct("SURFACE", "data/surface/tms/tms.gro")
      boxes[sim_label].add_mol(short, "data/mols/"+mol+"/"+mol+".gro", "fill", mass = mass[short], auto_dens=dens[short])
      boxes[sim_label].add_topol(pore_link+"pore.top", "master")
      boxes[sim_label].add_topol(pore_link+"grid.itp", "top")
      boxes[sim_label].add_topol(["data/surface/tms/tms.top", "data/surface/tms/tmsg.itp", "data/mols/"+mol+"/"+mol+".top"])
      boxes[sim_label].set_job(job)
      boxes[sim_label].set_param(param)
      boxes[sim_label].add_charge_si(1.28)


Create simulation objects
-------------------------

.. code-block:: python

  # Simulation series
  sim_series = ps.Simulate("sim/amorph_paper_res", [box for name, box in boxes.items()])
  sim_series.set_cluster(cluster)

  # Single simulation
  sim_single = ps.Simulate("output/single", box["1OL"])
  sim_single.set_cluster(cluster)


Generate folder structure
-------------------------

.. code-block:: python

  sim_series.generate()
  sim_single.generate()


Create benchmark
----------------

.. code-block:: python

  bench = ps.Benchmark(box["1OL"], 75, list(range(1,20+1)), "output/bench")
  bench.generate()


.. raw:: html

        </div>
      </div>
    </div>
  </div>
