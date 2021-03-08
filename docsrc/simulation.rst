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


Set individual job and parameter dictionaries
---------------------------------------------

.. code-block:: python

  job = {"min": {"file": "data/simulation/forhlr.sh", "nodes": 2, "np": 20, "wall": "24:00:00"},
         "nvt": {"file": "data/simulation/forhlr.sh", "nodes": 4, "np": 20, "wall": "24:00:00"},
         "run": {"file": "data/simulation/forhlr.sh", "maxh": 24, "nodes": 11, "np": 20, "runs": 15, "wall": "24:00:00"}}

  param = {"min": {"file": "data/simulation/pore_min.mdp"},
           "nvt": {"file": "data/simulation/pore_nvt.mdp", "param": {"NUMBEROFSTEPS": 2000000, "TEMPERATURE_VAL": 298}},
           "run": {"file": "data/simulation/pore_run.mdp", "param": {"NUMBEROFSTEPS": 20000000, "TEMPERATURE_VAL": 298}}}


Create box objects
------------------

.. code-block:: python

  box1 = ps.Box("box1")
  box1.add_box("data/simulation/pore.gro")
  box1.add_pore("data/simulation/pore.obj")
  box1.add_mol("EDC", "data/simulation/educt.gro", 10)
  box1.add_mol("PRD", "data/simulation/productmc.gro", 12)
  box1.add_mol("BEN", "data/simulation/benzene.gro", "fill", auto_dens=500)
  box1.add_topol("data/simulation/pore.top", "master")
  box1.add_topol("data/simulation/grid.itp", "top")
  box1.add_topol(["data/simulation/educt.top", "data/simulation/productmc.top", "data/simulation/benzene.top"])
  box1.add_topol(["data/simulation/tms.top", "data/simulation/tmsg.itp"])
  box1.set_job(job)
  box1.set_param(param)
  box1.add_charge_si(1.314730)

  box2 = ps.Box("box2", "bxx")
  box2.add_box("data/simulation/pore.gro")
  box2.add_pore("data/simulation/pore.obj")
  box2.add_mol("EDC", "data/simulation/educt.gro", 15)
  box2.add_mol("PRD", "data/simulation/productmc.gro", 12)
  box2.add_mol("BEN", "data/simulation/benzene.gro", "fill", auto_dens=500)
  box2.add_topol("data/simulation/pore.top", "master")
  box2.add_topol("data/simulation/grid.itp", "top")
  box2.add_topol(["data/simulation/educt.top", "data/simulation/productmc.top", "data/simulation/benzene.top"])
  box2.add_topol(["data/simulation/tms.top", "data/simulation/tmsg.itp"])
  box2.add_charge_si(1.314730)


Create simulation objects
-------------------------

.. code-block:: python

  sim1 = ps.Simulate("output/series", [box1, box2])  # Series
  sim2 = ps.Simulate("output/single", box1)  # Single


Generate folder structure
-------------------------

.. code-block:: python

  sim1.generate()
  sim2.generate()


Create Benchmark
----------------

.. code-block:: python

  bench = ps.Benchmark(box1,20,list(range(1,20+1)),"output/bench")
  bench.generate()


.. raw:: html

        </div>
      </div>
    </div>
  </div>
