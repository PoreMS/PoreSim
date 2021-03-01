################################################################################
# Simulate Class                                                               #
#                                                                              #
"""Create simulation boxes with all necessary files."""
################################################################################


import os
import yaml
import shutil

import poresim.utils as utils

from poresim.box import Box
from poresim.parameter import Parameter
from poresim.topology import Topology
from poresim.construct import Construct
from poresim.actuate import Actuate
from poresim.analyze import Analyze


class Simulate:
    """This class generates the simulation folder for a series of simulations.

    The output of this class can then be uploaded to a computing cluster.

    Parameters
    ----------
    link : string, optional
        Simulation folder link
    box : Box, list, optional
        Simulation box object, a series of those
    """
    def __init__(self, link="./simulation/", box=None):
        # Initialize
        self._link = link if link[-1] == "/" else link+"/"
        self._config = self._link+"/"+"config.yml"

        # Create simulation folder if not existent
        utils.mkdirp(self._link)

        # Create simulation config file if not existent
        if not os.path.isfile(self._config):
            # Get package directory
            package_dir = os.path.split(__file__)[0]+"/"

            # Job file
            file_in = package_dir+"templates/config.yml"
            file_out = self._config
            shutil.copy(file_in, file_out)

        # Load simulation config file
        with open(self._config, "r") as stream:
            self._sim_dict = yaml.load(stream, Loader=yaml.FullLoader)

            # Create box objects
            boxes = {}
            for i in self._sim_dict["box"]:
                temp = Box()
                temp.set_sim_dict(self._sim_dict["box"][i])
                boxes[i] = temp
            self._sim_dict["box"] = boxes

        # Add boxes to simulation dictionary
        if box is not None:
            boxes = box if isinstance(box, list) else [box]
            self._sim_dict["box"] = {i:boxes[i] for i in range(len(boxes))}
            self._update_config()


    ###################
    # Private Methods #
    ###################
    def _update_config(self):
        """Update simulation config file.
        """
        sim_dict = self._sim_dict.copy()

        sim_dict["box"] = {box_id:sim_dict["box"][box_id].get_sim_dict() for box_id in sim_dict["box"]}

        with open(self._config, "w") as file_out:
            file_out.write(yaml.dump(sim_dict))


    ##################
    # Public Methods #
    ##################
    def add_box(self, box):
        """Add simulation box to simulation dictionary.

        Parameters
        ----------
        box : Box, list
            Simulation box or a list of those
        """
        boxes = box if isinstance(box, list) else [box]
        for box in boxes:
            self._sim_dict["box"][len(boxes)] = box
        self._update_config()

    def generate(self):
        """Generate simulation folder.
        """
        # Run through boxes
        for box_id in self._sim_dict["box"]:
            box = self._sim_dict["box"][box_id]
            box_link = self._link+box.get_name()+"/" if len(self._sim_dict["box"])>1 else self._link

            # Create folder if multiple boxes exist
            utils.mkdirp(box_link)

            # Create parameter files
            params = self._sim_dict["param"] if box.get_param() is None else box.get_param()
            param = Parameter(box_link, params)
            param.generate_files()

            # Create topology files
            topol = Topology(box_link, box.get_topol())
            topol.generate_files()

            # Create structure files and shells
            construct = Construct(self._link, box_link, box.get_mols(), box.get_struct())
            construct.generate_files()

            # Create simulation shells
            job = self._sim_dict["job"] if box.get_job() is None else box.get_job()
            actuate = Actuate(self._link, box_link, self._sim_dict["cluster"], job, box.get_label(), box.get_struct())
            actuate.generate_files()

            # Create analysis shell file for automated filling
            for mol in box.get_mols():
                if box.get_mols()[mol][0]=="fill":
                    analyze = Analyze(self._link, box_link)
                    analyze.extract_mol("nvt")
                    # Copy automated filling
                    if box.get_mols()[mol][2] is not None:
                        utils.copy(os.path.split(__file__)[0]+"/templates/auto_dens.py", box_link+"ana/ana.py")
                        utils.replace(box_link+"ana/ana.py", "MOLSHORT", mol)
                        utils.replace(box_link+"ana/ana.py", "MOLLINK", "../_gro/"+box.get_struct()[mol].split("/")[-1])
                        utils.replace(box_link+"ana/ana.py", "TARGETDENS", str(box.get_mols()[mol][2]))
                        utils.replace(box_link+"ana/ana.py", "SUBMITCOMMAND", self._sim_dict["cluster"]["queuing"]["submit"]+" min.job")

            # End message
            print("Finished simulation folder - "+box.get_label()+" ...")


    ##################
    # Setter Methods #
    ##################
    def set_sim_dict(self, sim_dict):
        """Set the simulation dictionary.

        Parameters
        ----------
        sim_dict : dictionary
            Simulation dictionary
        """
        self._sim_dict = sim_dict
        self._update_config()

    def set_box(self, box):
        """Set the simulation box.

        Parameters
        ----------
        box : Box, list
            Simulation box object or a series of those
        """
        boxes = box if isinstance(box, list) else [box]
        self._sim_dict["box"] = {i:boxes[i] for i in range(len(boxes))}
        self._update_config()

    def set_job(self, job):
        """Set the dictionary containing all simulation information.

        Parameters
        ----------
        job : dictionary
            Information dictionary

        Examples
        --------
        Following dictionary structure may be used as an input

        .. code-block:: python

            {"min": {"file": "data/simulation/forhlr.sh", "nodes": 2, "np": 20, "wall": "24:00:00"},
             "nvt": {"file": "data/simulation/forhlr.sh", "nodes": 4, "np": 20, "wall": "24:00:00"},
             "run": {"file": "data/simulation/forhlr.sh", "maxh": 24, "nodes": 11, "np": 20, "runs": 15, "wall": "24:00:00"}}
        """
        self._sim_dict["job"] = job
        self._update_config()

    def set_param(self, param):
        """Set the dictionary containing all systems parameters.

        Parameters
        ----------
        param : dictionary
            Information dictionary

        Examples
        --------
        Following dictionary structure may be used as an input

        .. code-block:: python

            {"min": {"file": "data/simulation/pore_min.mdp"},
             "nvt": {"file": "data/simulation/pore_nvt.mdp", "param": {"NUMBEROFSTEPS": 2000000, "TEMPERATURE_VAL": 298}},
             "run": {"file": "data/simulation/pore_run.mdp", "param": {"NUMBEROFSTEPS": 20000000, "TEMPERATURE_VAL": 298}}}
        """
        self._sim_dict["param"] = param
        self._update_config()

    def set_cluster(self, cluster):
        """Set the dictionary containing all cluster preferences.

        Parameters
        ----------
        cluster : dictionary
            Information dictionary

        Examples
        --------
        Following dictionary structure may be used as an input

        .. code-block:: python

            {"address": "user_name@cluster",
             "directory": "/home/pores/simulation/",
             "queuing": {"add_np": False, "mpi": "$DO_PARALLEL", "shell": "forhlr.sh", "submit": "sbatch --partition multinode"}}
        """
        self._sim_dict["cluster"] = cluster
        self._update_config()


    ##################
    # Getter Methods #
    ##################
    def get_sim_dict(self):
        """Return the simulation dictionary.

        Returns
        -------
        sim_dict : dictionary
            Simulation dictionary
        """
        return self._sim_dict

    def get_box(self):
        """Return the simulation box dictionary.

        Returns
        -------
        box : dictionary
            Simulation box dictionary
        """
        return self._sim_dict["box"]
