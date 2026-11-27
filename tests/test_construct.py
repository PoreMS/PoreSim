"""Tests for the box construction shells and the GROMACS position files.

Every test builds its own simulation folder in a temporary directory and
drives :class:`poresim.construct.Construct` with a hand-written molecule
dictionary, so each placement branch (reservoir, pore, wall, box, areas,
fixed numbers and fills) is exercised on its own and its output checked.
"""

import os

import numpy as np
import pytest
import yaml

import poresim as ps
from poresim.construct import Construct

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
BENZENE = os.path.join(DATA, "benzene.gro")


def _mol(
    inp,
    *,
    auto_dens=None,
    mass=None,
    section="both",
    area=None,
    box=None,
    kwargs_gmx=None,
):
    """Build a molecule entry in the layout of Box.add_mol."""
    return [
        inp,
        12,
        auto_dens,
        mass,
        section,
        area or [],
        box or [],
        kwargs_gmx or {},
    ]


def _pore_yml(tmp_path, extra_shapes=None):
    """Write a pore description based on tests/data/pore.yml.

    ``extra_shapes`` maps shape ids to shape entries that are added next to the
    cylindrical shape_00 of the reference file.
    """
    with open(os.path.join(DATA, "pore.yml")) as f:
        props = yaml.safe_load(f)
    props.update(extra_shapes or {})
    path = tmp_path / "pore.yml"
    with open(path, "w") as f:
        yaml.dump(props, f)
    return str(path), props


def _slit_shape(props):
    """A slit shape along z with the dimensions of the reference cylinder."""
    cyl = props["shape_00"]
    return {
        "diameter": 3.0,
        "parameter": {
            "central": [0, 0, 1],
            "centroid": cyl["parameter"]["centroid"],
            "height": 3.0,
            "length": cyl["parameter"]["length"],
        },
        "shape": "SLIT",
    }


def _construct(tmp_path, mols, pore=None):
    """Run Construct.generate_files and return the folder and construct.sh."""
    link = str(tmp_path / "sim") + "/"
    os.makedirs(link)
    struct = dict.fromkeys(mols, BENZENE)
    struct["BOX"] = BENZENE
    if pore is not None:
        struct["PORE"] = pore
    Construct(link, link, mols, struct).generate_files()
    with open(link + "construct.sh") as f:
        return link, f.read()


def _positions(path):
    return np.loadtxt(path, ndmin=2)


def _read(path):
    with open(path) as f:
        return f.read()


#################
# Pore systems  #
#################
def test_pore_fixed_number_split_over_pore_and_reservoir(tmp_path):
    pore, _ = _pore_yml(tmp_path)
    _, construct = _construct(tmp_path, {"MOL": _mol(7, section="both")}, pore=pore)
    assert "# Fill Reservoir MOL" in construct
    assert "# Fill Pore MOL" in construct
    assert "-nmol 7" in construct
    # Pore index and topology update are written for pore systems
    assert "gmx_mpi make_ndx" in construct
    assert 'countMOL=$(($(grep -c "MOL" ./_gro/box.gro)/12))' in construct


def test_pore_fixed_number_split_over_the_shapes(tmp_path):
    """A fixed number in the pore is split evenly over the shapes along z."""
    shape_01 = _shape_along_x()
    shape_01["parameter"]["central"] = [0, 0, 1]
    pore, _ = _pore_yml(tmp_path, {"shape_01": shape_01})
    link, _ = _construct(tmp_path, {"MOL": _mol(7, section="pore")}, pore=pore)
    # 7 // 2 = 3 for the first shape, the remainder goes to the last one
    assert len(_positions(link + "_gro/position_shape_00_MOL.dat")) == 3
    assert len(_positions(link + "_gro/position_shape_01_MOL.dat")) == 4


def test_pore_fixed_number_skips_shapes_not_along_z(tmp_path):
    pore, _ = _pore_yml(tmp_path, {"shape_01": _shape_along_x()})
    link, construct = _construct(tmp_path, {"MOL": _mol(8, section="pore")}, pore=pore)
    assert os.path.isfile(link + "_gro/position_shape_00_MOL.dat")
    assert not os.path.isfile(link + "_gro/position_shape_01_MOL.dat")
    assert construct.count("-ip ./_gro/position_shape_00_MOL.dat") == 1
    assert "position_shape_01_MOL.dat" not in construct


def _shape_along_x():
    """A cylindrical shape along x, which the pore placement ignores."""
    with open(os.path.join(DATA, "pore.yml")) as f:
        shape = yaml.safe_load(f)["shape_00"]
    shape["parameter"]["central"] = [1, 0, 0]
    return shape


def test_pore_wall_layer(tmp_path):
    """section="wall" puts the molecules as two layers onto each slit wall."""
    _, props = _pore_yml(tmp_path)
    slit_x = _slit_shape(props)
    slit_x["parameter"] = {**slit_x["parameter"], "central": [1, 0, 0]}
    pore, props = _pore_yml(
        tmp_path, {"shape_00": _slit_shape(props), "shape_01": slit_x}
    )
    link, construct = _construct(tmp_path, {"MOL": _mol(10, section="wall")}, pore=pore)
    # Only slits along z get an insertion command
    assert "Filled shape_01 MOL" not in construct

    positions = _positions(link + "_gro/position_shape_00_MOL.dat")
    assert len(positions) == 10
    centroid = props["shape_00"]["parameter"]["centroid"]
    half_h = 3.0 / 2 * 0.75
    assert set(np.round(positions[:, 1], 6)) == {
        round(centroid[1] - half_h, 6),
        round(centroid[1] + half_h, 6),
    }
    assert "# Fill Pore Wall MOL" in construct
    assert "-nmol 10" in construct
    assert 'echo "Filled shape_00 MOL ..."' in construct


def test_pore_wall_fill_inserts_nothing_up_front(tmp_path):
    """A wall fill starts with zero molecules and leaves the rest to fill.sh."""
    _, props = _pore_yml(tmp_path)
    pore, _ = _pore_yml(tmp_path, {"shape_00": _slit_shape(props)})
    link, construct = _construct(
        tmp_path,
        {"MOL": _mol("fill", auto_dens=800, mass=78.11, section="wall")},
        pore=pore,
    )
    assert "# Fill Pore Wall MOL" in construct
    assert "-nmol 0" in construct
    assert os.path.isfile(link + "_fill/fill.sh")


@pytest.mark.parametrize("section", ["pore", "both"])
def test_pore_fill_by_density(tmp_path, section):
    # shape_01 lies along x and gets no molecules
    pore, props = _pore_yml(tmp_path, {"shape_01": _shape_along_x()})
    dens, mass = 800, 78.11
    link, construct = _construct(
        tmp_path,
        {"MOL": _mol("fill", auto_dens=dens, mass=mass, section=section)},
        pore=pore,
    )

    shape = props["shape_00"]
    factor = dens / mass / 10 * 6.022
    expected = int(
        factor * np.pi * shape["diameter"] ** 2 / 4 * shape["parameter"]["length"]
    )
    assert len(_positions(link + "_gro/position_shape_00_MOL.dat")) == expected
    assert f"-nmol {expected}" in construct
    assert os.path.isfile(link + "_gro/position_MOL.dat") == (section == "both")
    assert not os.path.isfile(link + "_gro/position_shape_01_MOL.dat")

    fill = _read(link + "_fill/fill.sh")
    assert "-nmol FILLDENS_MOL" in fill
    assert "mv ../_gro/index.ndx" in fill
    assert os.path.isfile(link + "_fill/fillBackup.sh")
    assert os.path.isfile(link + "_fill/sort.py")


def test_pore_fill_in_slit(tmp_path):
    """The slit volume is width x height x length instead of a cylinder."""
    _, props = _pore_yml(tmp_path)
    pore, props = _pore_yml(tmp_path, {"shape_00": _slit_shape(props)})
    dens, mass = 800, 78.11
    _, construct = _construct(
        tmp_path,
        {"MOL": _mol("fill", auto_dens=dens, mass=mass, section="pore")},
        pore=pore,
    )
    dims = props["system"]["dimensions"]
    slit = props["shape_00"]
    expected = int(
        dens
        / mass
        / 10
        * 6.022
        * dims[0]
        * slit["diameter"]
        * slit["parameter"]["length"]
    )
    assert f"-nmol {expected}" in construct
    assert f"-dr {0.9 * dims[0] / 2} {0.9 * slit['diameter'] / 2}" in construct


FILL_WITHOUT_DENSITY = pytest.mark.xfail(
    raises=TypeError,
    strict=True,
    reason=(
        "Open: Box.add_mol accepts a fill with a mass but without auto_dens, and "
        "fill.sh handles auto_dens=None, but Construct divides auto_dens by the "
        "mass to estimate the first fill and fails with a TypeError"
    ),
)


@FILL_WITHOUT_DENSITY
def test_pore_fill_without_density_refills_a_fixed_count(tmp_path):
    pore, _ = _pore_yml(tmp_path)
    link, _ = _construct(
        tmp_path, {"MOL": _mol("fill", mass=78.11, section="res")}, pore=pore
    )
    fill = _read(link + "_fill/fill.sh")
    assert "-nmol 0 " in fill
    # Without a target density there is no analysis to back up
    assert "# Backup Analysis" not in fill
    assert not os.path.isfile(link + "_fill/fillBackup.sh")


################
# Box systems  #
################
def test_box_fixed_number(tmp_path):
    _, construct = _construct(tmp_path, {"MOL": _mol(25, box=[4, 4, 4])})
    assert "# Fill Box" in construct
    assert "-nmol 25" in construct
    # No pore, so no index file
    assert "make_ndx" not in construct


def test_box_fixed_number_needs_box_dimensions(tmp_path, capsys):
    link = str(tmp_path) + "/"
    struct = {"MOL": BENZENE}
    mols = {"MOL": _mol(25, area=[[0, 2]], box=[4, 4, 4])}
    c = Construct(link, link, {"MOL": _mol(25)}, struct)
    os.makedirs(link + "_gro")
    c._pos_dat()
    assert "specify in add_mol for every molecule the box dimension" in (
        capsys.readouterr().out
    )

    c = Construct(link, link, mols, struct)
    c._pos_dat()
    assert len(_positions(link + "_gro/position_MOL_area0.dat")) == 25


def test_box_fixed_number_next_to_a_fill(tmp_path):
    """With a fill present, the fixed molecules get a position file as well."""
    link, _ = _construct(
        tmp_path,
        {
            "FIX": _mol(6, box=[4, 4, 4]),
            "MOL": _mol(
                "fill", auto_dens=800, mass=78.11, section="box", box=[4, 4, 4]
            ),
        },
    )
    positions = _positions(link + "_gro/position_FIX.dat")
    assert len(positions) == 6
    np.testing.assert_allclose(positions[0], [2, 2, 2])


def test_box_fill_by_density(tmp_path):
    dens, mass, box = 800, 78.11, [4, 4, 6]
    link, construct = _construct(
        tmp_path,
        {"MOL": _mol("fill", auto_dens=dens, mass=mass, section="box", box=box)},
    )
    expected = int(dens / mass / 10 * 6.022 * box[0] * box[1] * box[2])
    assert f"-nmol {expected}" in construct
    assert f"-dr {box[0] / 2} {box[1] / 2} {box[2] / 2}" in construct

    positions = _positions(link + "_gro/position_MOL.dat")
    assert len(positions) == int(
        dens / mass / 10 * 6.022 * box[0] * box[1] * box[2] * 0.8
    )
    np.testing.assert_allclose(positions[0], [2, 2, 3])

    fill = _read(link + "_fill/fill.sh")
    assert f"-dr {box[0] / 2} {box[1] / 2} {box[2] / 2}" in fill
    assert "index.ndx" not in fill


def test_box_fill_by_density_in_areas(tmp_path):
    dens, mass, box = 800, 78.11, [4, 4, 10]
    area = [[0, 3], [7, 10]]
    link, construct = _construct(
        tmp_path,
        {
            "MOL": _mol(
                "fill", auto_dens=dens, mass=mass, section="box", area=area, box=box
            )
        },
    )
    factor = dens / mass / 10 * 6.022
    for i, seg in enumerate(area):
        positions = _positions(link + f"_gro/position_MOL_area{i}.dat")
        assert len(positions) == int(factor * box[0] * box[1] * (seg[1] - seg[0]) * 0.8)
        assert positions[0][2] == pytest.approx((seg[0] + seg[1]) / 2)
        expected = int(factor * box[0] * box[1] * (seg[1] - seg[0]))
        assert f"-ip ./_gro/position_MOL_area{i}.dat -nmol {expected}" in construct

    fill = _read(link + "_fill/fill.sh")
    assert fill.count("-nmol FILLDENS_MOL") == 2


@FILL_WITHOUT_DENSITY
def test_box_fill_in_areas_without_density(tmp_path):
    link, _ = _construct(
        tmp_path,
        {"MOL": _mol("fill", mass=78.11, area=[[0, 3]], box=[4, 4, 10])},
    )
    assert "-nmol 10000" in _read(link + "_fill/fill.sh")


def test_gmx_kwargs_are_appended(tmp_path):
    _, construct = _construct(
        tmp_path,
        {"MOL": _mol(5, box=[4, 4, 4], kwargs_gmx={"-try": 50, "-scale": 0.4})},
    )
    assert "-nmol 5 -try 50 -scale 0.4 >> logging.log" in construct


def test_construct_in_box_subfolder(tmp_path):
    """A box in a subfolder of the simulation refers to it as ./<box>/."""
    sim = str(tmp_path) + "/"
    box = sim + "box1/"
    os.makedirs(box)
    struct = {"MOL": BENZENE, "BOX": BENZENE}
    Construct(sim, box, {"MOL": _mol(3, box=[4, 4, 4])}, struct).generate_files()
    construct = _read(sim + "construct.sh")
    assert "# Process ./box1/ #" in construct
    assert "-f ./box1/_gro/box.gro" in construct


############################
# Simulate for box systems #
############################
def _box_simulation_job():
    shell = os.path.join(DATA, "forhlr.sh")
    step = {"file": shell, "nodes": 1, "np": 4, "wall": "1:00:00"}
    return {"min": dict(step), "nvt": dict(step), "run": {**step, "maxh": 1, "runs": 2}}


def _box_simulation_param():
    param = {"NUMBEROFSTEPS": 10, "TEMPERATURE_VAL": 298}
    return {
        "min": {"file": os.path.join(DATA, "pore_min.mdp")},
        "nvt": {"file": os.path.join(DATA, "pore_nvt.mdp"), "param": dict(param)},
        "run": {"file": os.path.join(DATA, "pore_run.mdp"), "param": dict(param)},
    }


def _box_simulation(tmp_path, **mol_kwargs):
    box = ps.Box("box")
    box.add_box(os.path.join(DATA, "box.gro"))
    box.add_mol("BEN", BENZENE, **mol_kwargs)
    box.add_topol(os.path.join(DATA, "benzene.top"), "master")
    box.set_job(_box_simulation_job())
    box.set_param(_box_simulation_param())
    link = str(tmp_path / "sim")
    ps.Simulate(link, box).generate()
    return link + "/"


def test_simulate_box_fill_in_areas_renders_the_area_analysis(tmp_path):
    link = _box_simulation(
        tmp_path,
        inp="fill",
        auto_dens=800,
        mass=78.11,
        section="box",
        area=[[0, 2], [4, 6]],
        box=[3, 3, 6],
    )
    ana = _read(link + "ana/ana.py")
    assert "import poreana" in ana
    # The fill analysis of a box system evaluates every area on its own
    for low, up in ([0, 2], [4, 6]):
        assert f".index({low})" in ana
        assert f".index({up})" in ana
    assert 'position_{}_area{}.dat".format("BEN",0)' in ana
    assert 'position_{}_area{}.dat".format("BEN",1)' in ana
    # Fill systems are analysed after the NVT equilibration
    assert "-f ../nvt/nvt.trr" in _read(link + "ana/ana.sh")


def test_simulate_box_without_fill_is_analysed_after_the_run(tmp_path):
    link = _box_simulation(tmp_path, inp=20, section="box", box=[3, 3, 6])
    assert os.path.isfile(link + "ana/ana.py")
    assert "-f ../run/run.trr" in _read(link + "ana/ana.sh")


def test_simulate_box_fill_without_areas(tmp_path):
    link = _box_simulation(
        tmp_path, inp="fill", auto_dens=800, mass=78.11, section="box", box=[3, 3, 6]
    )
    ana = _read(link + "ana/ana.py")
    # Without areas the whole box is evaluated and refilled at its centre
    assert 'position_{}.dat".format("BEN")' in ana
    assert "str(3/2)" in ana


def test_simulate_reloads_an_existing_folder(tmp_path):
    link = _box_simulation(tmp_path, inp=20, section="box", box=[3, 3, 6])
    sim = ps.Simulate(link)
    boxes = sim.get_box()
    assert len(boxes) == 1
    assert boxes[0].get_mols()["BEN"][0] == 20


def test_simulate_pore_without_fill_extracts_nothing(tmp_path):
    box = ps.Box("box")
    box.add_box(os.path.join(DATA, "pore.gro"))
    box.add_pore(os.path.join(DATA, "pore.yml"))
    box.add_mol("BEN", BENZENE, 10)
    box.add_topol(os.path.join(DATA, "benzene.top"), "master")
    link = str(tmp_path / "sim") + "/"
    sim = ps.Simulate(link, box)
    sim.set_job(_box_simulation_job())
    sim.set_param(_box_simulation_param())
    sim.generate()
    assert os.path.isfile(link + "ana/ana.py")
    assert not os.path.isfile(link + "ana/ana.sh")
