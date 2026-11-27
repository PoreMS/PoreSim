import os
import tempfile

import numpy as np
import pytest

import poresim as ps
from poresim.construct import Construct
from poresim.topology import Topology


def test_utils():
    assert round(ps.utils.mumol_m2_to_mols(3, 100), 4) == 180.66
    assert round(ps.utils.mols_to_mumol_m2(180, 100), 4) == 2.989
    assert round(ps.utils.mmol_g_to_mumol_m2(0.072, 512), 2) == 0.14
    assert round(ps.utils.mmol_l_to_mols(30, 1000), 4) == 18.066
    assert round(ps.utils.mols_to_mmol_l(18, 1000), 4) == 29.8904


def test_box_api():
    """Test Box getter/setter round-trips."""
    box = ps.Box("mybox", "mylabel")
    assert box.get_name() == "mybox"
    assert box.get_label() == "mylabel"
    assert box.get_job() is None
    assert box.get_param() is None
    assert box.get_mols() == {}
    assert box.get_struct() == {}
    assert box.get_topol() == {"master": [], "top": [], "itp": []}

    box.set_name("newname")
    assert box.get_name() == "newname"

    box.set_label("newlabel")
    assert box.get_label() == "newlabel"

    job = {
        "run": {
            "file": "data/forhlr.sh",
            "nodes": 1,
            "np": 4,
            "wall": "1:00:00",
            "runs": 2,
        }
    }
    box.set_job(job)
    assert box.get_job() == job

    param = {"run": {"file": "data/pore_run.mdp"}}
    box.set_param(param)
    assert box.get_param() == param

    # sim_dict round-trip
    sd = box.get_sim_dict()
    box2 = ps.Box()
    box2.set_sim_dict(sd)
    assert box2.get_name() == "newname"
    assert box2.get_label() == "newlabel"


def test_box_add_mol_validation():
    """Test that add_mol rejects invalid inputs."""
    box = ps.Box()
    # float inp
    assert box.add_mol("X", "data/benzene.gro", 0.5) is None
    # bad num_atoms
    assert box.add_mol("X", "data/benzene.gro", 10, num_atoms="DOTA") is None
    # bad auto_dens
    assert box.add_mol("X", "data/benzene.gro", 10, auto_dens="heavy") is None
    # fill without mass
    assert box.add_mol("X", "data/benzene.gro", "fill", auto_dens=500) is None


def test_box_topol():
    """Test add_topol for all three types."""
    box = ps.Box()
    box.add_topol("data/pore.top", "master")
    box.add_topol(["data/educt.top", "data/benzene.top"])
    box.add_topol("data/grid.itp", "top")
    topol = box.get_topol()
    assert "data/pore.top" in topol["master"]
    assert "data/educt.top" in topol["itp"]
    assert "data/benzene.top" in topol["itp"]
    assert "data/grid.itp" in topol["top"]


def test_construct_write_positions():
    """Test that _write_positions writes correct x y z lines."""
    struct = {"MOL": "data/benzene.gro"}
    mols = {"MOL": [10, 12, None, None, "both", [], [], {}]}
    c = Construct("output/", "output/", mols, struct)

    positions = np.array([[1.5, 2.5, 3.5], [4.0, 5.0, 6.0], [0.1, 0.2, 0.3]])
    with tempfile.NamedTemporaryFile(mode="r", suffix=".dat", delete=False) as f:
        path = f.name

    try:
        c._write_positions(path, positions)
        with open(path) as f:
            lines = f.readlines()
        assert len(lines) == 3
        parts0 = lines[0].strip().split()
        assert len(parts0) == 3
        assert float(parts0[0]) == pytest.approx(1.5)
        assert float(parts0[1]) == pytest.approx(2.5)
        assert float(parts0[2]) == pytest.approx(3.5)
        parts2 = lines[2].strip().split()
        assert float(parts2[0]) == pytest.approx(0.1)
    finally:
        os.unlink(path)


def test_construct_gmx_insert_cmd():
    """Test that _gmx_insert_cmd builds correct GROMACS commands."""
    struct = {"MOL": "data/benzene.gro"}
    mols = {"MOL": [10, 12, None, None, "both", [], [], {}]}
    c = Construct("output/", "output/", mols, struct)

    # Basic command
    cmd = c._gmx_insert_cmd("../_gro/", "box.gro", "MOL", 100)
    assert "gmx_mpi insert-molecules" in cmd
    assert "-f ../_gro/box.gro" in cmd
    assert "-o ../_gro/box.gro" in cmd
    assert "-ci ../_gro/benzene.gro" in cmd
    assert "-nmol 100" in cmd
    assert ">> logging.log 2>&1" in cmd
    assert "-dr" not in cmd
    assert "-ip" not in cmd

    # With dr
    cmd_dr = c._gmx_insert_cmd("../_gro/", "box.gro", "MOL", 50, dr=[1.0, 2.0, 3.0])
    assert "-dr 1.0 2.0 3.0" in cmd_dr

    # With pos_file
    cmd_ip = c._gmx_insert_cmd("../_gro/", "box.gro", "MOL", 50, pos_file="pos.dat")
    assert "-ip pos.dat" in cmd_ip

    # With both
    cmd_both = c._gmx_insert_cmd(
        "../_gro/", "box.gro", "MOL", 50, dr=[0.5, 0.5, 5.0], pos_file="pos.dat"
    )
    assert "-dr 0.5 0.5 5.0" in cmd_both
    assert "-ip pos.dat" in cmd_both

    # With kwargs_gmx
    struct2 = {"MOL": "data/benzene.gro"}
    mols2 = {"MOL": [10, 12, None, None, "both", [], [], {"-try": 500, "-scale": 0.5}]}
    c2 = Construct("output/", "output/", mols2, struct2)
    cmd_kw = c2._gmx_insert_cmd("../_gro/", "box.gro", "MOL", 10)
    assert "-try 500" in cmd_kw
    assert "-scale 0.5" in cmd_kw


def test_box_add_mol_with_given_atom_number():
    box = ps.Box()
    box.add_mol("X", "data/benzene.gro", 10, num_atoms=12)
    assert box.get_mols()["X"][1] == 12


def _write(path, text):
    with open(path, "w") as f:
        f.write(text)
    return str(path)


def test_topology_itp_without_charge(tmp_path):
    """Without a silicon charge the CHARGESI placeholder stays in the itp."""
    top = _write(
        tmp_path / "mol.top",
        "[ defaults ]\n1 2 yes 0.5 0.8333\n\n"
        "[ moleculetype ]\nMOL 3\n\n"
        "[ atoms ]\n1 SI 1 MOL SI1 1 CHARGESI 28.0855\n\n"
        "[ system ]\nMOL\n",
    )
    Topology(str(tmp_path) + "/", {"itp": [top]}).generate_files()
    itp = (tmp_path / "_top" / "mol.itp").read_text()
    assert "[ defaults ]" not in itp
    assert "[ system ]" not in itp
    assert "[ moleculetype ]" in itp
    assert "CHARGESI" in itp


def test_topology_itp_with_charge(tmp_path):
    top = _write(
        tmp_path / "grid.top", "[ atoms ]\n1 SI 1 MOL SI1 1 CHARGESI 28.0855\n\n"
    )
    Topology(str(tmp_path) + "/", {"itp": [top], "charge": 1.28}).generate_files()
    assert "1.280000" in (tmp_path / "_top" / "grid.itp").read_text()


@pytest.mark.xfail(
    strict=True,
    reason=(
        "Open: Topology._read keeps only the first block of a section that "
        "occurs twice, so e.g. the improper dihedrals that acpype writes as a "
        "second [ dihedrals ] block are dropped from the itp"
    ),
)
def test_topology_keeps_repeated_sections(tmp_path):
    top = _write(
        tmp_path / "mol.top",
        "[ dihedrals ] ; propers\n1 2 3 4 9\n\n"
        "[ dihedrals ] ; impropers\n1 2 3 5 4\n\n",
    )
    Topology(str(tmp_path) + "/", {"itp": [top]}).generate_files()
    itp = (tmp_path / "_top" / "mol.itp").read_text()
    assert "1 2 3 4 9" in itp
    assert "1 2 3 5 4" in itp


def test_benchmark_without_npt(tmp_path):
    """Without an npt step the benchmark only adjusts min and nvt."""
    shell = "data/forhlr.sh"
    step = {"file": shell, "nodes": 1, "np": 4, "wall": "1:00:00"}
    job = {"min": dict(step), "nvt": dict(step), "run": {**step, "maxh": 1, "runs": 2}}
    param = {
        "min": {"file": "data/pore_min.mdp"},
        "nvt": {
            "file": "data/pore_nvt.mdp",
            "param": {"NUMBEROFSTEPS": 10, "TEMPERATURE_VAL": 298},
        },
        "run": {
            "file": "data/pore_run.mdp",
            "param": {"NUMBEROFSTEPS": 10, "TEMPERATURE_VAL": 298},
        },
    }
    box = ps.Box("box")
    box.add_box("data/box.gro")
    box.add_mol("BEN", "data/benzene.gro", 5, box=[3, 3, 3])
    box.add_topol("data/benzene.top", "master")
    box.set_job(job)
    box.set_param(param)

    link = str(tmp_path / "bench")
    bench = ps.Benchmark(box, [4], [1, 2], link)
    bench.set_job(job)
    bench.set_param(param)
    bench.generate()

    shell_text = (tmp_path / "bench" / "benchmark.sh").read_text()
    assert "cp -r X/_gro 01/_gro" in shell_text
    assert "cp -r X/_gro 02/_gro" in shell_text
    assert not (tmp_path / "bench" / "01" / "npt").exists()
