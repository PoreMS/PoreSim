import os
import shutil

import pytest
import poresim as ps


@pytest.fixture(autouse=True, scope="module")
def clean_output():
    """Create and clean the output directory before integration tests run."""
    folder = "output"
    ps.utils.mkdirp(folder)
    ps.utils.mkdirp(folder + "/temp")
    open(folder + "/temp.txt", "a").close()
    for filename in os.listdir(folder):
        file_path = os.path.join(folder, filename)
        if os.path.isfile(file_path) or os.path.islink(file_path):
            os.unlink(file_path)
        elif os.path.isdir(file_path):
            shutil.rmtree(file_path)
    yield


@pytest.fixture(scope="module")
def job():
    return {
        "min": {"file": "data/forhlr.sh", "nodes": 2, "np": 20, "wall": "24:00:00"},
        "nvt": {"file": "data/forhlr.sh", "nodes": 4, "np": 20, "wall": "24:00:00"},
        "run": {"file": "data/forhlr.sh", "maxh": 24, "nodes": 11, "np": 20, "runs": 15, "wall": "24:00:00"},
    }


@pytest.fixture(scope="module")
def param():
    return {
        "min": {"file": "data/pore_min.mdp"},
        "nvt": {"file": "data/pore_nvt.mdp", "param": {"NUMBEROFSTEPS": 2000000, "TEMPERATURE_VAL": 298}},
        "run": {"file": "data/pore_run.mdp", "param": {"NUMBEROFSTEPS": 20000000, "TEMPERATURE_VAL": 298}},
    }


def _check_position_file(path):
    """Assert a position file exists and every line has exactly 3 floats."""
    assert os.path.isfile(path), f"Position file missing: {path}"
    with open(path) as f:
        lines = f.readlines()
    assert len(lines) > 0, f"Empty position file: {path}"
    for line in lines:
        parts = line.strip().split()
        assert len(parts) == 3, f"Line does not have 3 columns in {path}: {line!r}"
        for p in parts:
            float(p)  # raises ValueError if not a number


def test_utils_file_io():
    file_link = "output/test/test.txt"
    ps.utils.mkdirp("output/test")
    with open(file_link, "w") as f:
        f.write("TEST")
    ps.utils.copy(file_link, file_link + "t")
    ps.utils.replace(file_link + "t", "TEST", "DOTA")
    with open(file_link + "t", "r") as f:
        for line in f:
            assert line == "DOTA\n"


def test_box(job, param):
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

    assert box.add_mol("EDC", "data/educt.gro", 0.1) is None
    assert box.add_mol("EDC", "data/educt.gro", 10, num_atoms="DOTA") is None
    assert box.add_mol("EDC", "data/educt.gro", 10, auto_dens="DOTA") is None


def test_fill_box(job, param):
    box = ps.Box("box")
    box.add_box("data/box.gro")
    box.add_mol("BEN", "data/benzene.gro", "fill", auto_dens=500)
    box.add_topol("data/benzene.top", "master")
    box.set_job(job)
    box.set_param(param)


def test_sim(job, param):
    cluster = {
        "address": "user_name@cluster",
        "directory": "/home/pores/simulation/",
        "queuing": {
            "add_np": False,
            "mpi": "$DO_PARALLEL",
            "shell": "forhlr.sh",
            "submit": "sbatch --partition multinode",
        },
    }

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

    sim2 = ps.Simulate("output/single", box1)
    sim3 = ps.Simulate("output/single", box1)
    sim3.add_box(box2)
    assert len(sim3.get_box()) == 2, "add_box should increment box count"
    sim3.set_sim_dict(sim2.get_sim_dict())
    sim3.set_box(list(sim2.get_box().values()))
    sim3.set_cluster(cluster)

    sim2.generate()

    base = "output/single"
    # Required shell scripts exist
    assert os.path.isfile(f"{base}/construct.sh")
    assert os.path.isfile(f"{base}/equilibrate.sh")
    assert os.path.isfile(f"{base}/simulate.sh")

    # construct.sh contains GROMACS insert-molecules and topology update
    with open(f"{base}/construct.sh") as f:
        construct = f.read()
    assert "gmx_mpi insert-molecules" in construct
    assert "count" in construct  # topology update loop

    # Job shells for each step exist and contain GROMACS commands
    for step in ("min", "nvt", "run", "run0"):
        path = f"{base}/{step.replace('0', '')}/{step}.job" if step != "run0" else f"{base}/run/run0.job"
        assert os.path.isfile(path), f"Missing job file: {path}"
        with open(path) as f:
            content = f.read()
        assert "gmx_mpi" in content

    # Position files exist and have correct format
    gro = f"{base}/_gro"
    _check_position_file(f"{gro}/position_BEN.dat")
    _check_position_file(f"{gro}/position_shape_00_BEN.dat")
    _check_position_file(f"{gro}/position_EDC.dat")
    _check_position_file(f"{gro}/position_shape_00_EDC.dat")
    _check_position_file(f"{gro}/position_PRD.dat")
    _check_position_file(f"{gro}/position_shape_00_PRD.dat")

    # fill.sh exists and contains refill logic
    assert os.path.isfile(f"{base}/_fill/fill.sh")
    with open(f"{base}/_fill/fill.sh") as f:
        fill = f.read()
    assert "gmx_mpi insert-molecules" in fill
    assert "fill_num" in fill

    # ana/ana.py exists and was rendered from template
    assert os.path.isfile(f"{base}/ana/ana.py")
    with open(f"{base}/ana/ana.py") as f:
        ana = f.read()
    assert "import poreana" in ana


def test_bench(job, param):
    npt_job = dict(job)
    npt_job["npt"] = {"file": "data/forhlr.sh", "nodes": 4, "np": 20, "wall": "24:00:00"}
    npt_param = dict(param)
    npt_param["npt"] = {"file": "data/pore_nvt.mdp", "param": {"NUMBEROFSTEPS": 2000000, "TEMPERATURE_VAL": 298}}

    box = ps.Box("box")
    box.add_box("data/pore.gro")
    box.add_pore("data/pore.yml")
    box.add_mol("EDC", "data/educt.gro", 10, section="pore")
    box.add_mol("PRD", "data/productmc.gro", 12, section="res")
    box.add_mol("BEN", "data/benzene.gro", "fill", auto_dens=500, mass=78.11)
    box.add_topol("data/pore.top", "master")
    box.add_topol("data/grid.itp", "top")
    box.add_topol(["data/educt.top", "data/productmc.top", "data/benzene.top"])
    box.add_topol(["data/tms.top", "data/tmsg.itp"])
    box.add_struct("GRO", "data/benzene.gro")
    box.set_job(npt_job)
    box.set_param(npt_param)
    box.add_charge_si(1.28)

    bench1 = ps.Benchmark(box, 20, list(range(21)), "output/bench1")
    bench1.set_job(npt_job)
    bench1.set_param(npt_param)
    bench1.generate()

    bench2 = ps.Benchmark(box, 20, list(range(21)), "output/bench2", iterator="np")
    bench2.set_job(npt_job)
    bench2.set_param(npt_param)
    bench2.generate()

    # Verify benchmark folders were created with expected structure
    assert os.path.isdir("output/bench1")
    assert os.path.isfile("output/bench1/construct.sh")
    assert os.path.isfile("output/bench1/benchmark.sh")
    # Position files present in the pore box
    _check_position_file("output/bench1/X/_gro/position_BEN.dat")
    _check_position_file("output/bench1/X/_gro/position_shape_00_BEN.dat")


def test_2phase(job, param):
    pores = ps.Box("353_2phase")
    pores.set_label("353_2phase")
    pores.add_box("data/box_2phase.gro")
    pores.add_mol("CAT", "data/2phase/catalyst.gro", inp=10, area=[[0, 5], [15, 20]], box=[8, 8, 20], kwargs_gmx={"-try": 1000, "-scale": 0.47})
    pores.add_mol("EDC", "data/2phase/reactant.gro", inp=10, area=[[5, 15]], box=[8, 8, 20], kwargs_gmx={"-try": 1000, "-scale": 0.47})
    pores.add_mol("IM", "data/2phase/bmi.gro", inp=1600, area=[[0, 5], [15, 20]], box=[8, 8, 20], kwargs_gmx={"-try": 1000, "-scale": 0.47})
    pores.add_mol("BF4", "data/2phase/bf4.gro", inp=1620, area=[[0, 5], [15, 20]], box=[8, 8, 20], kwargs_gmx={"-try": 1000, "-scale": 0.47})
    pores.add_mol("HEP", "data/2phase/1-heptane.gro", inp=2420, area=[[5, 15]], box=[8, 8, 20], kwargs_gmx={"-try": 1000, "-scale": 0.47})
