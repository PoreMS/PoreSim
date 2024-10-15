################################################################################
# Benchmark Class                                                             #
#                                                                              #
"""All necessary function for running the simulation."""
################################################################################


import poresim.utils as utils

from poresim.simulate import Simulate
from poresim.parameter import Parameter
from poresim.topology import Topology
from poresim.construct import Construct
from poresim.actuate import Actuate


class Benchmark(Simulate):
    """This class extends the simulate class
    :class:`poresim.simulate.Simulate` for benchmarking. It creates a
    benchmark in order to figure out the best possible number of node and
    processor combination to simulate with.

    Depending on the iterator, different aspects will be benchmarked. Currently
    following inputs are possible

    * **nodes** - varying number of nodes at constant CPU number
    * **np** - varying number of CPUs at constant node number

    Parameters
    ----------
    box : Box
        Simulation box object
    np : integer, list
        Number of processors to try
    nodes : integer, list
        Number of nodes to try
    link : string, optional
        Simulation folder link
    iterator : string, optional
        Aspect to be benchmarked
    """
    def __init__(self, box, np, nodes, link="./benchmark", iterator="nodes"):
        # Call super class
        super(Benchmark, self).__init__(link, box)

        # Initialize
        self._box = box
        self._np = np if isinstance(np, list) else [np]
        self._nodes = nodes if isinstance(nodes, list) else [nodes]
        self._iterator = iterator


    ##################
    # Public Methods #
    ##################
    def generate(self):
        """Generate simulation folder.
        """
        # Set parameter dictionary
        params = self._sim_dict["param"]
        params["nvt"]["param"]["NUMBEROFSTEPS"] = 5000
        if "npt" in params:
            params["npt"]["param"]["NUMBEROFSTEPS"] = 5000
        if "run" in params:
            del params["run"]

        # Remove run from job dictionary
        if "run" in self._sim_dict["job"]:
            del self._sim_dict["job"]["run"]

        # Set box specifications
        self._box.set_name("X")
        self._box.set_label("X")
        self._box.set_param(params)

        # Create parameter files
        box_link = self._link+self._box.get_name()+"/"

        # Create folder if multiple boxes exist
        utils.mkdirp(box_link)

        # Create parameter files
        param = Parameter(box_link, self._box.get_param())
        param.generate_files()

        # Create topology files
        topol = Topology(box_link, self._box.get_topol())
        topol.generate_files()

        # Create structure files and shells
        construct = Construct(self._link, box_link, self._box.get_mols(), self._box.get_struct())
        construct.generate_files()

        # Create simulation shells
        actuate = Actuate(self._link, box_link, self._sim_dict["cluster"], self._sim_dict["job"], self._box.get_label(), self._box.get_struct())
        actuate.generate_files()

        # Process iterator input
        if self._iterator == "nodes":
            iterator = self._nodes
        elif self._iterator == "np":
            iterator = self._np

        # Create copy shell after minimization
        with open(self._link+"benchmark.sh", "a") as file_out:
            # Create subfolders for different simulations
            job = self._sim_dict["job"]
            for it in iterator:
                # Set nodes and np
                n = it if iterator == self._np else self._np[0]
                node = it if iterator == self._nodes else self._nodes[0]

                # Set job dictionary
                job["min"]["np"] = n
                job["min"]["nodes"] = node
                job["min"]["wall"] = "01:00:00"

                job["nvt"]["np"] = n
                job["nvt"]["nodes"] = node
                job["nvt"]["wall"] = "01:00:00"

                if "npt" in params:
                    job["npt"]["np"] = n
                    job["npt"]["nodes"] = node
                    job["npt"]["wall"] = "01:00:00"

                # Set subfolder
                self._sim_dict["job"] = job
                self._box.set_name(str(it).zfill(2))
                self._box.set_label("b_"+str(it).zfill(2))
                box_link = self._link+self._box.get_name()+"/"

                # Create parameter files
                param = Parameter(box_link, self._box.get_param())
                param.generate_files()

                # Create simulation shells
                actuate = Actuate(self._link, box_link, self._sim_dict["cluster"], self._sim_dict["job"], self._box.get_label(), self._box.get_struct())
                actuate.generate_files()

                # Write copy to shell
                file_out.write("cp -r X/_gro "+self._box.get_name()+"/_gro\n")
                file_out.write("cp -r X/_top "+self._box.get_name()+"/_top\n")
                file_out.write("cp X/min/* "+self._box.get_name()+"/min/\n\n")

        # Set equilibration master shell to nvt
        utils.replace(self._link+"equilibrate.sh", "min", "nvt")
