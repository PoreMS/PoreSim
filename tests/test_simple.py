import os
import sys

import shutil
import unittest

import poresim as ps


class UserModelCase(unittest.TestCase):
    #################
    # Remove Output #
    #################
    @classmethod
    def setUpClass(self):
        if os.path.isdir("tests"):
            os.chdir("tests")

        folder = 'output'
        ps.utils.mkdirp(folder)
        ps.utils.mkdirp(folder+"/temp")
        open(folder+"/temp.txt", 'a').close()

        for filename in os.listdir(folder):
            file_path = os.path.join(folder, filename)
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)

    #########
    # Utils #
    #########
    def test_utils(self):
        file_link = "output/test/test.txt"

        ps.utils.mkdirp("output/test")

        with open(file_link, "w") as file_out:
            file_out.write("TEST")
        ps.utils.copy(file_link, file_link+"t")
        ps.utils.replace(file_link+"t", "TEST", "DOTA")
        with open(file_link+"t", "r") as file_in:
            for line in file_in:
                self.assertEqual(line, "DOTA\n")

        self.assertEqual(round(ps.utils.mumol_m2_to_mols(3, 100), 4), 180.66)
        self.assertEqual(round(ps.utils.mols_to_mumol_m2(180, 100), 4), 2.989)
        self.assertEqual(round(ps.utils.mmol_g_to_mumol_m2(0.072, 512), 2), 0.14)
        self.assertEqual(round(ps.utils.mmol_l_to_mols(30, 1000), 4), 18.066)
        self.assertEqual(round(ps.utils.mols_to_mmol_l(18, 1000), 4), 29.8904)


    ############
    # Simulate #
    ############
    def test_box(self):
        job = {"min": {"file": "data/forhlr.sh", "nodes": 2, "np": 20, "wall": "24:00:00"},
               "nvt": {"file": "data/forhlr.sh", "nodes": 4, "np": 20, "wall": "24:00:00"},
               "run": {"file": "data/forhlr.sh", "maxh": 24, "nodes": 11, "np": 20, "runs": 15, "wall": "24:00:00"}}

        param = {"min": {"file": "data/pore_min.mdp"},
                 "nvt": {"file": "data/pore_nvt.mdp", "param": {"NUMBEROFSTEPS": 2000000, "TEMPERATURE_VAL": 298}},
                 "run": {"file": "data/pore_run.mdp", "param": {"NUMBEROFSTEPS": 20000000, "TEMPERATURE_VAL": 298}}}

        # Boxes
        box = ps.Box("box")
        box.add_box("data/pore.gro")
        box.add_pore("data/pore.yml")
        box.add_mol("EDC", "data/educt.gro", 10)
        box.add_mol("PRD", "data/productmc.gro", 12)
        box.add_mol("BEN", "data/benzene.gro", "fill", auto_dens=500)
        box.add_topol("data/pore.top", "master")
        box.add_topol("data/grid.itp", "top")
        box.add_topol(["data/educt.top", "data/productmc.top", "data/benzene.top"])
        box.add_topol(["data/tms.top", "data/tmsg.itp"])
        box.add_struct("GRO", "data/benzene.gro")
        box.set_job(job)
        box.set_param(param)
        box.add_charge_si(1.28)

        print()
        self.assertEqual(box.add_mol("EDC", "data/educt.gro", 0.1), None)
        self.assertEqual(box.add_mol("EDC", "data/educt.gro", 10, num_atoms="DOTA"), None)
        self.assertEqual(box.add_mol("EDC", "data/educt.gro", 10, auto_dens="DOTA"), None)


    def test_fill_box(self):
        job = {"min": {"file": "data/forhlr.sh", "nodes": 2, "np": 20, "wall": "24:00:00"},
               "nvt": {"file": "data/forhlr.sh", "nodes": 4, "np": 20, "wall": "24:00:00"},
               "run": {"file": "data/forhlr.sh", "maxh": 24, "nodes": 11, "np": 20, "runs": 15, "wall": "24:00:00"}}

        param = {"min": {"file": "data/pore_min.mdp"},
                 "nvt": {"file": "data/pore_nvt.mdp", "param": {"NUMBEROFSTEPS": 2000000, "TEMPERATURE_VAL": 298}},
                 "run": {"file": "data/pore_run.mdp", "param": {"NUMBEROFSTEPS": 20000000, "TEMPERATURE_VAL": 298}}}

        # Boxes
        box = ps.Box("box")
        box.add_box("data/box.gro")
        box.add_mol("BEN", "data/benzene.gro", "fill", auto_dens=500)
        box.add_topol("data/benzene.top", "master")
        box.set_job(job)
        box.set_param(param)

        print()


    ############
    # Simulate #
    ############
    def test_sim(self):
        job = {"min": {"file": "data/forhlr.sh", "nodes": 2, "np": 20, "wall": "24:00:00"},
               "nvt": {"file": "data/forhlr.sh", "nodes": 4, "np": 20, "wall": "24:00:00"},
               "run": {"file": "data/forhlr.sh", "maxh": 24, "nodes": 11, "np": 20, "runs": 15, "wall": "24:00:00"}}

        param = {"min": {"file": "data/pore_min.mdp"},
                 "nvt": {"file": "data/pore_nvt.mdp", "param": {"NUMBEROFSTEPS": 2000000, "TEMPERATURE_VAL": 298}},
                 "run": {"file": "data/pore_run.mdp", "param": {"NUMBEROFSTEPS": 20000000, "TEMPERATURE_VAL": 298}}}

        cluster = {"address": "user_name@cluster",
                   "directory": "/home/pores/simulation/",
                   "queuing": {"add_np": False, "mpi": "$DO_PARALLEL", "shell": "forhlr.sh", "submit": "sbatch --partition multinode"}}

        box1 = ps.Box("box1")
        box1.add_box("data/pore.gro")
        box1.add_pore("data/pore.yml")
        box1.add_mol("EDC", "data/educt.gro", 10)
        box1.add_mol("PRD", "data/productmc.gro", 12)
        box1.add_mol("BEN", "data/benzene.gro", "fill", auto_dens=500, mass=78.11)
        box1.add_topol("data/pore.top", "master")
        box1.add_topol("data/grid.itp", "top")
        box1.add_topol(["data/educt.top", "data/productmc.top", "data/benzene.top"])
        box1.add_topol(["data/tms.top", "data/tmsg.itp"])
        box1.add_struct("GENERATE", "data/benzene.gro")
        box1.add_struct("PLUMED", "data/benzene.gro")
        box1.set_job(job)
        box1.set_param(param)
        box1.add_charge_si(1.314730)

        box2 = ps.Box("box2", "bxx")
        box2.add_box("data/pore.gro")
        box2.add_pore("data/pore.yml")
        box2.add_mol("EDC", "data/educt.gro", 15)
        box2.add_mol("PRD", "data/productmc.gro", 12)
        box2.add_mol("BEN", "data/benzene.gro", "fill", auto_dens=500, mass=78.11)
        box2.add_topol("data/pore.top", "master")
        box2.add_topol("data/grid.itp", "top")
        box2.add_topol(["data/educt.top", "data/productmc.top", "data/benzene.top"])
        box2.add_topol(["data/tms.top", "data/tmsg.itp"])
        box2.set_job(job)
        box2.set_param(param)
        box2.add_charge_si(1.314730)

        #sim1 = ps.Simulate("output/series", [box1, box2])  # Series
        sim2 = ps.Simulate("output/single", box1)  # Single
        sim3 = ps.Simulate("output/single", box1)

        sim3.add_box(box2)
        sim3.set_sim_dict(sim2.get_sim_dict())
        sim3.set_box(list(sim2.get_box().values()))
        sim3.set_cluster(cluster)

        print()
        #sim1.generate()
        print()
        sim2.generate()


    #############
    # Benchmark #
    #############
    def test_bench(self):
        job = {"min": {"file": "data/forhlr.sh", "nodes": 2, "np": 20, "wall": "24:00:00"},
               "nvt": {"file": "data/forhlr.sh", "nodes": 4, "np": 20, "wall": "24:00:00"},
               "npt": {"file": "data/forhlr.sh", "nodes": 4, "np": 20, "wall": "24:00:00"},
               "run": {"file": "data/forhlr.sh", "maxh": 24, "nodes": 11, "np": 20, "runs": 15, "wall": "24:00:00"}}

        param = {"min": {"file": "data/pore_min.mdp"},
                 "nvt": {"file": "data/pore_nvt.mdp", "param": {"NUMBEROFSTEPS": 2000000, "TEMPERATURE_VAL": 298}},
                 "npt": {"file": "data/pore_nvt.mdp", "param": {"NUMBEROFSTEPS": 2000000, "TEMPERATURE_VAL": 298}},
                 "run": {"file": "data/pore_run.mdp", "param": {"NUMBEROFSTEPS": 20000000, "TEMPERATURE_VAL": 298}}}

        box = ps.Box("box")
        box.add_box("data/pore.gro")
        box.add_pore("data/pore.yml")
        box.add_mol("EDC", "data/educt.gro", 10)
        box.add_mol("PRD", "data/productmc.gro", 12)
        box.add_mol("BEN", "data/benzene.gro", "fill", auto_dens=500, mass=78.11)
        box.add_topol("data/pore.top", "master")
        box.add_topol("data/grid.itp", "top")
        box.add_topol(["data/educt.top", "data/productmc.top", "data/benzene.top"])
        box.add_topol(["data/tms.top", "data/tmsg.itp"])
        box.add_struct("GRO", "data/benzene.gro")
        box.set_job(job)
        box.set_param(param)
        box.add_charge_si(1.28)

        bench1 = ps.Benchmark(box, 20, list(range(21)), "output/bench1")
        bench1.set_job(job)
        bench1.set_param(param)
        bench1.generate()

        bench2 = ps.Benchmark(box, 20, list(range(21)), "output/bench2", iterator="np")
        bench2.set_job(job)
        bench2.set_param(param)
        bench2.generate()


if __name__ == '__main__':
    unittest.main(verbosity=2)
