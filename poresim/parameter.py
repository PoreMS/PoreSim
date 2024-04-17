################################################################################
# Parameter Class                                                              #
#                                                                              #
"""Process parameter files."""
################################################################################


import poresim.utils as utils


class Parameter:
    """This class generates parameter files by replacing the given parameter in
    prepared files.

    Parameters
    ----------
    link : string
        Simulation box link
    param : dictionary
        Parameter dictionary
    """
    def __init__(self, link, param):
        # Initialize
        self._link = link+"_mdp/"
        self._param = param

        # Create folder
        utils.mkdirp(self._link)

    ##################
    # Public Methods #
    ##################
    def generate_files(self):
        """Create a parameter file for simulation based on the given parameter
        dictionary of the simulation box.

        The information for each simulation step should contain a file link, and
        optionally a list of parameters, to be replaced in the file. Parameter
        have the syntax of replacing string followed by the value.

        Examples
        --------
        >>> TEMPERATURE_VAL: 298
        """
        # Run through mdp files
        for step in self._param:
            file_name = step+".mdp"
            utils.copy(self._param[step]["file"], self._link+file_name)

            if "param" in self._param[step]:
                for inp in self._param[step]["param"]:
                    utils.replace(self._link+file_name, inp, str(self._param[step]["param"][inp]))
