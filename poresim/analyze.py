################################################################################
# Analyse Class                                                                #
#                                                                              #
"""All necessary function for creating analysis shells."""
################################################################################


import poresim.utils as utils


class Analyze:
    """This class creates shell-files for generating shell file for analysing
    a simulation.

    Parameters
    ----------
    sim_link : string
        Simulation master folder link
    box_link : string
        Simulation box folder link
    """
    def __init__(self, sim_link, box_link, cluster):
        # Initialize
        self._sim_link = sim_link
        self._box_path = box_link
        self._cluster = cluster
        self._box_link = "./" if sim_link == box_link else "./" + \
            box_link.split("/")[-2]+"/"


    ##################
    # Public Methods #
    ##################
    def extract_mol(self, step):
        """This function extracts molecules of a type from a trajectory for
        further analysis using MolDyn.

        Parameters
        ----------
        step : string
            Step to analyze (**nvt**, **npt**, **run**)
        """
        # Set folder names
        folder_step = "../"+step+"/"

        # Set file names
        file_trr = step+".trr"
        file_tpr = step+".tpr"

        # Open file
        utils.mkdirp(self._box_path+"ana")
        with open(self._box_path+"ana/ana.sh", "w") as file_out:
            # Check if backup folder is given
            file_out.write("# Set Todos\n")
            file_out.write("echo \"Load gromacs ...\"; exit;\n")
            file_out.write("echo \"Set MOLECULEINDEX ...\"; exit;\n")
            file_out.write("\n")

            # Extract molecule
            file_out.write("# Extract molecules from trajectory\n")
            file_out.write("declare -A mols\n\n")

            file_out.write("mols[]=""\n\n")

            file_out.write("for key in \"${!mols[@]}\"; do\n")

            out_string = "gmx_mpi trjconv "
            out_string += "-f "+folder_step+file_trr+" "
            out_string += "-s "+folder_step+file_tpr+" "
            out_string += "-o traj_${mols[$key]}.xtc "
            out_string += "-pbc mol "
            out_string += "<<EOF\n"
            out_string += "$key\n"
            out_string += "EOF\n\n"
            file_out.write(out_string)
            
            file_out.write("done\n\n")


            file_out.write("echo \"System "+self._box_link+ " - Finished Extraction ...\"\n\n")

            file_out.write(self._cluster["queuing"]["submit"] + " ana.job")

            
