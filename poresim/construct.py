################################################################################
# Construct Class                                                              #
#                                                                              #
"""All necessary function for creating the finished simulation box."""
################################################################################


import os

import numpy as np
import poresim.utils as utils


class Construct:
    """This class creates shell-files for generating and filling the simulation
    box using GROMACS.

    Parameters
    ----------
    sim_link : string
        Simulation master folder link
    box_link : string
        Simulation box folder link
    mols : dictionary
        Molecules dictionary
    struct : dictionary
        Structure dictionary
    """
    def __init__(self, sim_link, box_link, mols, struct):
        self._sim_link = sim_link
        self._box_path = box_link
        self._box_link = "./" if sim_link == box_link else "./" + box_link.split("/")[-2] + "/"
        self._mols = mols
        self._struct = struct
        if "PORE" in struct:
            self._pore_props = utils.load(struct["PORE"])


    ###################
    # Private Methods #
    ###################
    def _write_positions(self, path, positions):
        """Write an (N, 3) numpy array of xyz coordinates to a Gromacs position file.

        Parameters
        ----------
        path : string
            Output file path
        positions : numpy.ndarray
            Array of shape (N, 3) with x, y, z coordinates
        """
        with open(path, "w", encoding="utf-8") as f:
            for row in positions:
                f.write(f"{row[0]} {row[1]} {row[2]}\n")

    def _gmx_insert_cmd(self, folder_gro, file_box, mol, nmol, dr=None, pos_file=None):
        """Build a gmx_mpi insert-molecules command string.

        Parameters
        ----------
        folder_gro : string
            Path to the gro folder
        file_box : string
            Box structure filename
        mol : string
            Molecule short name
        nmol : int or string
            Number of molecules to insert
        dr : list, optional
            Displacement radii [dx, dy, dz]
        pos_file : string, optional
            Position file path for the -ip flag

        Returns
        -------
        string
            Complete GROMACS command with log redirection
        """
        mol_file = self._struct[mol].split("/")[-1]
        parts = [
            "gmx_mpi insert-molecules",
            f"-f {folder_gro}{file_box}",
            f"-o {folder_gro}{file_box}",
            f"-ci {folder_gro}{mol_file}",
        ]
        if dr is not None:
            parts.append(f"-dr {dr[0]} {dr[1]} {dr[2]}")
        if pos_file is not None:
            parts.append(f"-ip {pos_file}")
        parts.append(f"-nmol {nmol}")
        for key, value in self._mols[mol][-1].items():
            parts.append(f"{key} {value}")
        return " ".join(parts) + " >> logging.log 2>&1\n"

    def _topol_index(self, file_out, path):
        """Update topology with molecule counts and create an index file for pore systems.

        Parameters
        ----------
        file_out : File
            Open file object to write into
        path : string
            Simulation root path (e.g. "./" or "./box/")
        """
        folder_gro = path + "_gro/"
        folder_top = path + "_top/"
        file_box = "box.gro"
        file_top = "topol.top"
        file_ndx = "index.ndx"

        file_out.write("# Update Topology\n")
        for mol in self._mols:
            file_out.write(f"count{mol}=$(($(grep -c \"{mol}\" {folder_gro}{file_box})/{self._mols[mol][1]}))\n")
            file_out.write(f"echo \"{mol} \"$count{mol} >> {folder_top}{file_top}\n")
        file_out.write(f"echo \"System {self._box_link} - Updated topology ...\"\n\n")

        if "PORE" in self._struct:
            file_out.write("# Create Index\n")
            cmd = (
                f"gmx_mpi make_ndx "
                f"-f {folder_gro}{file_box} "
                f"-o {folder_gro}{file_ndx}"
                f" >> logging.log 2>&1 <<EOF\n"
                "0 & a SI1 OM1\n"
                "q\n"
                "EOF\n"
            )
            file_out.write(cmd)
            file_out.write(f"echo \"System {self._box_link} - Created pore index file ...\"\n")

    def _pos_dat(self):
        """Create position files for Gromacs molecule insertion.

        Uses numpy for efficient generation of repeated coordinate arrays.
        Handles pore systems (reservoir/pore/wall sections) and box systems
        (fill by density or fixed number, with or without area constraints).
        """
        for mol in self._mols:
            inp, _, auto_dens, mass, section, area, box_dim, _ = self._mols[mol]
            is_fill = inp == "fill"
            gro = f"{self._box_path}_gro/"

            # Pore system: density-based fill (reservoir and/or pore regions)
            if is_fill and not box_dim and section in ("res", "both", "pore") and "PORE" in self._struct:
                res = self._pore_props["system"]["reservoir"]
                factor = auto_dens / mass / 10 * 6.022

                if section in ("res", "both"):
                    dims = self._pore_props["system"]["dimensions"]
                    num = int(factor * dims[0] * dims[1] * res * 2 * 0.5)
                    cx, cy = dims[0] / 2, dims[1] / 2
                    top = np.tile([cx, cy, res / 2], (num // 2, 1))
                    bot = np.tile([cx, cy, dims[2] - res / 2], (num // 2, 1))
                    self._write_positions(f"{gro}position_{mol}.dat", np.vstack([top, bot]))

                for pore_id in self._pore_props:
                    if pore_id[:5] != "shape":
                        continue
                    params = self._pore_props[pore_id]["parameter"]
                    if params["central"] != [0, 0, 1] or section not in ("pore", "both"):
                        continue
                    centroid = params["centroid"]
                    num_pore = int(factor * np.pi * self._pore_props[pore_id]["diameter"] ** 2 / 4 * params["length"])
                    positions = np.tile([centroid[0], centroid[1], res + centroid[2]], (num_pore, 1))
                    self._write_positions(f"{gro}position_{pore_id}_{mol}.dat", positions)

            # Pore system: fixed number of molecules
            elif not is_fill and not area and not box_dim and "PORE" in self._struct:
                j = sum(1 for pid in self._pore_props if pid[:5] == "shape")

                if section == "both":
                    num = [inp // (j + 2)] * (j + 2)
                elif section == "res":
                    num = [inp // 2] * 2
                elif section in ("pore", "box"):
                    num = [inp // j] * j
                else:
                    num = []

                if num and sum(num) != inp:
                    num[-1] += abs(inp - sum(num))

                if section in ("pore", "both"):
                    for pore_id, j in zip(self._pore_props, range(j)):
                        if pore_id[:5] != "shape":
                            continue
                        if self._pore_props[pore_id]["parameter"]["central"] != [0, 0, 1]:
                            continue
                        centroid = self._pore_props[pore_id]["parameter"]["centroid"]
                        res = self._pore_props["system"]["reservoir"]
                        positions = np.tile([centroid[0], centroid[1], res + centroid[2]], (num[j], 1))
                        self._write_positions(f"{gro}position_{pore_id}_{mol}.dat", positions)

                if section in ("res", "both"):
                    dims = self._pore_props["system"]["dimensions"]
                    res = self._pore_props["system"]["reservoir"]
                    top = np.tile([dims[0] / 2, dims[1] / 2, res / 2], (num[0], 1))
                    bot = np.tile([dims[0] / 2, dims[1] / 2, dims[2] - res / 2], (num[-1], 1))
                    self._write_positions(f"{gro}position_{mol}.dat", np.vstack([top, bot]))

                elif section == "wall":
                    for pore_id, j in zip(self._pore_props, range(j)):
                        if pore_id[:5] != "shape":
                            continue
                        p = self._pore_props[pore_id]["parameter"]
                        res = self._pore_props["system"]["reservoir"]
                        cx, cy = p["centroid"][0], p["centroid"][1]
                        cz = res + p["centroid"][2]
                        half_h = p["height"] / 2 * 0.75
                        n_half = inp // 2
                        top = np.tile([cx, cy - half_h, cz], (n_half, 1))
                        bot = np.tile([cx, cy + half_h, cz], (n_half, 1))
                        self._write_positions(f"{gro}position_{pore_id}_{mol}.dat", np.vstack([top, bot]))

            # Box system: density-based fill with area sections
            elif is_fill and area and section != "wall":
                factor = auto_dens / mass / 10 * 6.022
                for i, seg in enumerate(area):
                    num = int(factor * box_dim[0] * box_dim[1] * (seg[1] - seg[0]) * 0.8)
                    positions = np.tile([box_dim[0] / 2, box_dim[1] / 2, (seg[1] + seg[0]) / 2], (num, 1))
                    self._write_positions(f"{gro}position_{mol}_area{i}.dat", positions)

            # Box system: density-based fill without area
            elif is_fill and not area and section != "wall":
                factor = auto_dens / mass / 10 * 6.022
                num = int(factor * box_dim[0] * box_dim[1] * box_dim[2] * 0.8)
                positions = np.tile([box_dim[0] / 2, box_dim[1] / 2, box_dim[2] / 2], (num, 1))
                self._write_positions(f"{gro}position_{mol}.dat", positions)

            # Box system: fixed number without area
            elif not is_fill and not area:
                if not box_dim:
                    print("If you fill a system with one molecule specify in add_mol for every molecule the box dimension")
                    return
                positions = np.tile([box_dim[0] / 2, box_dim[1] / 2, box_dim[2] / 2], (inp, 1))
                self._write_positions(f"{gro}position_{mol}.dat", positions)

            # Box system: fixed number with area sections
            elif not is_fill and area:
                n_each = inp // len(area)
                for i, seg in enumerate(area):
                    positions = np.tile([box_dim[0] / 2, box_dim[1] / 2, (seg[1] + seg[0]) / 2], (n_each, 1))
                    self._write_positions(f"{gro}position_{mol}_area{i}.dat", positions)

    def _structure(self):
        """Create a shell file for constructing and filling the simulation box using GROMACS.
        Also updates the topology and (for pore systems) creates the required index file.
        """
        folder_gro = self._box_link + "_gro/"
        folder_fill = self._box_link + "_fill/"
        file_box = "box.gro"

        with open(self._sim_link + "construct.sh", "a", encoding="utf-8") as file_out:
            label_width = 12 + len(self._box_link)
            file_out.write("#" * label_width + "\n")
            file_out.write(f"# Process {self._box_link} #\n")
            file_out.write("#" * label_width + "\n")
            file_out.write("echo \"Load gromacs ...\"; exit;\n")
            if "fill" in [self._mols[mol][0] for mol in self._mols]:
                file_out.write("echo \"Set ions names in sort script if necessary ...\"; exit;\n")

            for mol in self._mols:
                inp, _, auto_dens, mass, section, area, box_dim, _ = self._mols[mol]
                is_fill = inp == "fill"
                file_out.write(f"\n####### {mol} #########\n")

                if "PORE" in self._struct:
                    dims = self._pore_props["system"]["dimensions"]
                    res = self._pore_props["system"]["reservoir"]
                    num = int(auto_dens / mass / 10 * 6.022 * dims[0] ** 2 * res * 2) if (is_fill and mass) else 0

                    if section in ("res", "both"):
                        file_out.write(f"# Fill Reservoir {mol}\n")
                        nmol = str(int(inp)) if not is_fill else str(num)
                        dr = [dims[0] / 2, dims[1] / 2, res / 2]
                        file_out.write(self._gmx_insert_cmd(folder_gro, file_box, mol, nmol, dr=dr, pos_file=f"{folder_gro}position_{mol}.dat"))
                        file_out.write(f"echo \"Filled reservoir with {mol} ...\"\n\n")

                    if section in ("pore", "both"):
                        file_out.write(f"# Fill Pore {mol}\n")
                        for pore_id in self._pore_props:
                            if pore_id[:5] != "shape":
                                continue
                            p = self._pore_props[pore_id]
                            if p["parameter"]["central"] != [0, 0, 1]:
                                continue
                            d = p["diameter"]
                            length = p["parameter"]["length"]
                            if is_fill and mass:
                                if p["shape"] == "SLIT":
                                    num_pore = int(auto_dens / mass / 10 * 6.022 * dims[0] * d * length)
                                else:
                                    num_pore = int(auto_dens / mass / 10 * 6.022 * np.pi * d ** 2 / 4 * length)
                            else:
                                num_pore = 0
                            if section != "res":
                                nmol = str(int(inp)) if not is_fill else str(num_pore)
                                if p["shape"] == "SLIT":
                                    dr = [0.90 * dims[0] / 2, 0.90 * d / 2, 0.9 * length / 2]
                                else:
                                    r_eff = 0.50 * np.sqrt(0.9 * d ** 2) / 2
                                    dr = [r_eff, r_eff, 0.9 * length / 2]
                                file_out.write(self._gmx_insert_cmd(folder_gro, file_box, mol, nmol, dr=dr, pos_file=f"{folder_gro}position_{pore_id}_{mol}.dat"))
                        file_out.write(f"echo \"Filled pore with {mol} ...\"\n\n")

                    elif section == "wall":
                        file_out.write(f"# Fill Pore Wall {mol}\n")
                        for pore_id in self._pore_props:
                            if pore_id[:5] != "shape":
                                continue
                            p = self._pore_props[pore_id]
                            if p["parameter"]["central"] != [0, 0, 1]:
                                continue
                            nmol = str(int(inp)) if not is_fill else "0"
                            dr = [0.9 * dims[0] / 2, p["diameter"] / 2 * 0.1, p["parameter"]["length"] / 2]
                            file_out.write(self._gmx_insert_cmd(folder_gro, file_box, mol, nmol, dr=dr, pos_file=f"{folder_gro}position_{pore_id}_{mol}.dat"))
                            file_out.write(f"echo \"Filled {pore_id} {mol} ...\"\n\n")

                else:
                    file_out.write("# Fill Box\n")
                    if area:
                        for i, seg in enumerate(area):
                            num = int(auto_dens / mass / 10 * 6.022 * box_dim[0] * box_dim[1] * (seg[1] - seg[0])) if is_fill else inp
                            nmol = str(int(inp)) if not is_fill else str(num)
                            dr = [box_dim[0] / 2, box_dim[1] / 2, (seg[1] - seg[0]) / 2]
                            file_out.write(self._gmx_insert_cmd(folder_gro, file_box, mol, nmol, dr=dr, pos_file=f"{folder_gro}position_{mol}_area{i}.dat"))
                    else:
                        num = int(auto_dens / mass / 10 * 6.022 * box_dim[0] * box_dim[1] * box_dim[2]) if is_fill else inp
                        if is_fill:
                            dr = [box_dim[0] / 2, box_dim[1] / 2, box_dim[2] / 2]
                            pos_file = f"{folder_gro}position_{mol}.dat"
                        else:
                            dr, pos_file = None, None
                        file_out.write(self._gmx_insert_cmd(folder_gro, file_box, mol, str(int(num)), dr=dr, pos_file=pos_file))

            file_out.write(f"python {folder_fill}sort.py {folder_gro}\n")
            file_out.write(f"echo \"System {self._box_link} - Filled simulation box ...\"\n\n")
            self._topol_index(file_out, self._box_link)
            file_out.write(f"rm {folder_gro}*#\n")
            file_out.write("rm logging.log\n\n")

    def _fill(self):
        """Create a shell file for continuously refilling a simulation box.

        The last equilibration steps are moved to a backup folder, the latest
        structure is restored, and the box is refilled via gmx insert-molecules.
        """
        sim_min = "min"
        sim_nvt = "nvt"

        folder_fill = "./"
        folder_gro = "../_gro/"
        folder_top = "../_top/"
        folder_ana = "../ana/"
        folder_min = f"../{sim_min}/"
        folder_nvt = f"../{sim_nvt}/"

        file_box = "box.gro"
        file_top = "topol.top"
        file_t_b = "topolBackup.top"
        file_ndx = "index.ndx"

        has_auto_dens = not all(self._mols[mol][2] is None for mol in self._mols)

        with open(self._box_path + "_fill/fill.sh", "w", encoding="utf-8") as file_out:
            file_out.write("# Create Todos\n")
            file_out.write("echo \"Load gromacs ...\"; exit;\n")
            if has_auto_dens:
                file_out.write("echo \"Load gromacs in Backup...\"; exit;\n")
            file_out.write("echo \"Set ions names in sort script if necessary ...\"; exit;\n\n")

            file_out.write("# Set folder number\n")
            file_out.write("fill_num=1\n\n")

            file_out.write("# Backup Simulation\n")
            file_out.write(f"mkdir {folder_fill}$fill_num\n")
            file_out.write(f"mv {folder_gro}{file_box} {folder_fill}$fill_num\n")
            file_out.write(f"mv {folder_top}{file_top} {folder_fill}$fill_num\n")
            if "PORE" in self._struct:
                file_out.write(f"mv {folder_gro}{file_ndx} {folder_fill}$fill_num\n")
            file_out.write(f"cp {folder_nvt}{sim_nvt}.gro {folder_gro}{file_box}\n")
            file_out.write(f"cp {folder_top}{file_t_b} {folder_top}{file_top}\n")
            file_out.write(f"mv {folder_min} {folder_fill}$fill_num\n")
            file_out.write(f"mv {folder_nvt} {folder_fill}$fill_num\n")
            file_out.write(f"mkdir {folder_min}\n")
            file_out.write(f"mkdir {folder_nvt}\n")
            file_out.write(f"cp {folder_fill}$fill_num/{sim_min}/{sim_min}.job {folder_min}\n")
            file_out.write(f"cp {folder_fill}$fill_num/{sim_nvt}/{sim_nvt}.job {folder_nvt}\n")
            file_out.write(f"echo \"System {self._box_link} - Backed up equilibration ...\"\n\n")

            if has_auto_dens:
                file_out.write("# Backup Analysis\n")
                file_out.write(f"mv {folder_ana} {folder_fill}$fill_num\n")
                file_out.write(f"mkdir {folder_ana}\n")
                file_out.write(f"cp {folder_fill}$fill_num/ana/ana.* {folder_ana}\n")
                file_out.write(f"echo \"System {self._box_link} - Backed up analysis ...\"\n\n")

            file_out.write("# Refill Box\n")
            for mol in self._mols:
                inp, _, auto_dens, mass, section, area, box_dim, _ = self._mols[mol]
                if inp != "fill":
                    continue

                mol_file = self._struct[mol].split("/")[-1]
                gmx_prefix = (
                    f"gmx_mpi insert-molecules "
                    f"-f {folder_gro}{file_box} "
                    f"-o {folder_gro}{file_box} "
                    f"-ci {folder_gro}{mol_file} "
                    f"-try 1000 -scale 0.47 "
                )

                if area:
                    for i, seg in enumerate(area):
                        nmol = "10000" if auto_dens is None else f"FILLDENS_{mol}"
                        dr = [box_dim[0] / 2, box_dim[1] / 2, (seg[1] - seg[0]) / 2]
                        file_out.write(
                            gmx_prefix +
                            f"-dr {dr[0]} {dr[1]} {dr[2]} "
                            f"-ip {folder_gro}position_{mol}_area{i}.dat "
                            f"-nmol {nmol} >> logging.log 2>&1\n"
                        )
                else:
                    nmol = "0" if auto_dens is None else f"FILLDENS_{mol}"
                    if "PORE" in self._struct:
                        dims = self._pore_props["system"]["dimensions"]
                        res = self._pore_props["system"]["reservoir"]
                        dr = [dims[0] / 2, dims[1] / 2, res / 2]
                    else:
                        dr = [box_dim[0] / 2, box_dim[1] / 2, box_dim[2] / 2]
                    file_out.write(
                        gmx_prefix +
                        f"-dr {dr[0]} {dr[1]} {dr[2]} "
                        f"-ip {folder_gro}position_{mol}.dat "
                        f"-nmol {nmol} >> logging.log 2>&1\n"
                    )

            file_out.write(f"python sort.py {folder_gro}\n")
            file_out.write(f"echo \"System {self._box_link} - Refilled simulation box ...\"\n\n")
            self._topol_index(file_out, "../")
            file_out.write("\n# Remove logs\n")
            file_out.write(f"rm {folder_gro}*#\n")
            file_out.write("rm logging.log\n\n")

            file_out.write("# Step fill folder number\n")
            file_out.write("cp fill.sh temp.sh\n")
            file_out.write("sed -i \"s/fill_num=$fill_num/fill_num=$((fill_num+1))/\" temp.sh\n")
            if has_auto_dens:
                file_out.write("sed -i \"s/fill_num=$fill_num/fill_num=$((fill_num+1))/\" fillBackup.sh\n")
            file_out.write("mv temp.sh fill.sh\n")


    ##################
    # Public Methods #
    ##################
    def generate_files(self):
        """Generate structure files and shells."""
        self._structure()

        utils.mkdirp(self._box_path + "_gro")

        for mol in self._struct:
            file_link = self._struct[mol]
            if mol == "BOX":
                utils.copy(file_link, self._box_path + "_gro/box.gro")
            elif mol == "GENERATE":
                utils.copy(file_link, self._box_path + "_gro/generate.sh")
            elif mol == "PLUMED":
                utils.copy(file_link, self._box_path + "_gro/plumed.dat")
            else:
                utils.copy(file_link, self._box_path + "_gro/" + file_link.split("/")[-1])

        if "PORE" in self._struct:
            sections = [self._mols[mol][4] for mol in self._mols]
            if ("wall" in sections or "pore" in sections) and "box" not in sections:
                self._pos_dat()

        for mol in self._mols:
            if self._mols[mol][5]:
                self._pos_dat()

        if "fill" in [self._mols[mol][0] for mol in self._mols]:
            utils.mkdirp(self._box_path + "_fill")
            self._pos_dat()
            self._fill()
            if not all(self._mols[mol][2] is None for mol in self._mols):
                utils.copy(self._box_path + "_fill/fill.sh", self._box_path + "_fill/fillBackup.sh")
            utils.copy(os.path.split(__file__)[0] + "/templates/sort.py", self._box_path + "_fill/sort.py")
