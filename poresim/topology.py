################################################################################
# Topology Class                                                               #
#                                                                              #
"""Process topology files."""
################################################################################


import poresim.utils as utils


class Topology:
    """This class reads and extract all information from topology files of the
    molecules specified in the simulation box. The aim was creating an
    independent class for editing and updating the topologies if necessary.

    Further editing functions should be easy to add.

    Parameters
    ----------
    link : string
        Simulation box link
    topol : dictionary
        Topology dictionary
    """
    def __init__(self, link, topol):
        # Initialize
        self._link = link+"_top/"
        self._topol = topol

        # Create folder
        utils.mkdirp(self._link)


    ###################
    # Private Methods #
    ###################
    def _read(self, link):
        """This function reads a topology file and saves all information and
        data in a dictionary, accessible with the section’s header.

        Parameters
        ----------
        link : string
            Topology-file link

        Returns
        -------
        data : dictionary
            Dictionary containing all data of the topology file
        """
        # Initialize
        data = {}

        # Read data
        is_new = False
        block = ""
        with open(link, "r") as file_in:
            for line in file_in:
                dat = line.split()

                if is_new:
                    if len(dat) > 0:
                        data[block].append(line)
                    else:
                        is_new = False

                elif len(dat) > 0 and dat[0] == "[":
                    block = dat[1]

                    if not dat[1] in data.keys():
                        is_new = True
                        data[block] = []

        return data

    def _itp(self, link, charge):
        """Create an itp file from a GROMACS topology contain everything except

        * defaults
        * atomtypes
        * system
        * molecules

        Parameters
        ----------
        link : string
            Topology file link
        charge : None, float
            Silicon charge for pore grid molecules
        """
        # Initialize
        data = self._read(link)
        dont = ["defaults", "atomtypes", "system", "molecules"]

        itp_link = self._link+link.split("/")[-1].split(".")[0]+".itp"

        # Create itp
        with open(itp_link, "w") as file_out:
            for param in data:
                if not param in dont:
                    file_out.write("[ "+param+" ]\n")
                    for line in data[param]:
                        file_out.write(line)
                    file_out.write("\n")

        # Replace charges in itp files
        if charge is not None:
            utils.replace(itp_link, "CHARGESI", "%8.6f" % charge)


    ##################
    # Public Methods #
    ##################
    def generate_files(self):
        """Generate topology files.
        """
        # Get silicon charge if existent
        charge = self._topol["charge"] if "charge" in self._topol else None

        # Run through topologies
        for top_type in self._topol:
            if not top_type=="charge":
                for file_link in self._topol[top_type]:
                    if top_type=="master":
                        utils.copy(file_link, self._link+"topol.top")
                        utils.copy(file_link, self._link+"topolBackup.top")
                    elif top_type=="top":
                        utils.copy(file_link, self._link+file_link.split("/")[-1])
                    elif top_type=="itp":
                        self._itp(file_link, charge)
