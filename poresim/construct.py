################################################################################
# Construct Class                                                              #
#                                                                              #
"""All necessary function for creating the finished simulation box."""
################################################################################


import os

import poresim.utils as utils


class Construct:
    """This class creates shell-files for generating the and filling the
    simulation box using GROMACS.

    Parameters
    ----------
    sim_link : string
        Simulation master folder link
    box_link : string
        Simulation box folder link
    struct : dictionary
        Structure dictionary
    """
    def __init__(self, sim_link, box_link, mols, struct):
        # Initialize
        self._sim_link = sim_link
        self._box_path = box_link
        self._box_link = "./" if sim_link == box_link else "./"+box_link.split("/")[-2]+"/"
        self._mols = mols
        self._struct = struct


    ###################
    # Private Methods #
    ###################
    def _topol_index(self, file_out, path):
        """Helper function for updating the topology and creating an index file,
        if the simulated system is a pore.

        Parameters
        ----------
        file_out : File
            File object to write in
        path : String
            Simulation root path
        """
        # Set folder names
        folder_gro = path+"_gro/"
        folder_top = path+"_top/"

        # Set file names
        file_box = "box.gro"
        file_top = "topol.top"
        file_ndx = "index.ndx"

        # Add number of residues to topology
        file_out.write("# Update Topology\n")
        for mol in self._mols:
            file_out.write("count"+mol+"=$(($(grep -c \""+mol+"\" "+folder_gro+file_box+")/"+str(self._mols[mol][1])+"))\n")
            file_out.write("echo \""+mol+" \"$count"+mol+" >> "+folder_top+file_top+"\n")

        file_out.write("echo \"System "+self._box_link+" - Updated topology ...\"\n\n")

        # Create index file
        if "PORE" in self._struct:
            file_out.write("# Create Index\n")
            out_string = "gmx_mpi make_ndx "
            out_string += "-f "+folder_gro+file_box+" "
            out_string += "-o "+folder_gro+file_ndx
            out_string += " >> logging.log 2>&1 <<EOF\n"
            out_string += "0 & a SI1 OM1\n"
            out_string += "q\n"
            out_string += "EOF\n"
            file_out.write(out_string)

            file_out.write("echo \"System "+self._box_link+" - Created pore index file ...\"\n")

    def _structure(self):
        """Create a shell file for constructing and filling the simulation
        box using GROMACS. Additionally, the master topology file is updated
        with the number of added molecules. In case a pore simulation is
        intended the needed index-file containing a new group for SI and OM is
        generated.
        """
        # Set folder names
        folder_gro = self._box_link+"_gro/"
        folder_fill = self._box_link+"_fill/"

        # Set file names
        file_box = "box.gro"

        # Open file
        with open(self._sim_link+"construct.sh", "a") as file_out:
            # Create box label
            file_out.write("#"*(12+len(self._box_link))+"\n")
            file_out.write("# Process "+self._box_link+" #\n")
            file_out.write("#"*(12+len(self._box_link))+"\n")
            file_out.write("echo \"Load gromacs ...\"; exit;\n")
            if "fill" in [self._mols[mol][0] for mol in self._mols]:
                file_out.write("echo \"Set ions names in sort script if necessary ...\"; exit;\n")

            # Fill box
            file_out.write("# Fill Box\n")
            for mol in self._mols:
                out_string = "gmx_mpi insert-molecules "
                out_string += "-f "+folder_gro+file_box+" "
                out_string += "-o "+folder_gro+file_box+" "
                out_string += "-ci "+folder_gro+self._struct[mol].split("/")[-1]+" "
                out_string += "-nmol "+str(int(self._mols[mol][0])) if not self._mols[mol][0]=="fill" else "-nmol "+str(10000)
                out_string += " >> logging.log 2>&1\n"
                if self._mols[mol][0]=="fill":
                    out_string += "python "+folder_fill+"empty_grid.py "+folder_gro+" "+mol+" "+str(self._mols[mol][1])+"\n"
                file_out.write(out_string)

            if "fill" in [self._mols[mol][0] for mol in self._mols]:
                file_out.write("python "+folder_fill+"sort.py "+folder_gro+"\n")

            file_out.write("echo \"System "+self._box_link+" - Filled simulation box ...\"\n\n")

            # Update topology and create index
            self._topol_index(file_out, self._box_link)

            # Remove log and backup
            file_out.write("rm "+folder_gro+"*#\n")
            file_out.write("rm logging.log\n\n")

    def _fill(self):
        """This function creates a shell file for continuously filling a pore
        simulation manually. The nvt and min equilibration of the last step are
        moved to a backup folder and the last structure is moved to the
        structure folder for the new construction and simulation.
        """
        # Set simulation folder descriptors
        sim_ana = "ana"
        sim_min = "min"
        sim_nvt = "nvt"

        # Set folder names
        folder_fill = "./"
        folder_gro = "../_gro/"
        folder_top = "../_top/"
        folder_ana = "../"+sim_ana+"/"
        folder_min = "../"+sim_min+"/"
        folder_nvt = "../"+sim_nvt+"/"

        # Set file names
        file_box = "box.gro"
        file_top = "topol.top"
        file_t_b = "topolBackup.top"
        file_ndx = "index.ndx"

        # Open file
        with open(self._box_path+"_fill/fill.sh", "w") as file_out:
            # Check if backup folder is given
            file_out.write("# Create Todos\n")
            file_out.write("echo \"Load gromacs ...\"; exit;\n")
            if not all(x is None for x in [self._mols[mol][2] for mol in self._mols]):
                file_out.write("echo \"Load gromacs in Backup...\"; exit;\n")
            file_out.write("echo \"Set ions names in sort script if necessary ...\"; exit;\n")
            file_out.write("\n")

            # Add folder number
            file_out.write("# Set folder number\n")
            file_out.write("fill_num=1\n\n")

            # Backup simulation
            file_out.write("# Backup Simulation\n")
            out_string = "mkdir "+folder_fill+"$fill_num\n"
            out_string += "mv "+folder_gro+file_box+" "+folder_fill+"$fill_num\n"
            out_string += "mv "+folder_top+file_top+" "+folder_fill+"$fill_num\n"
            out_string += "mv "+folder_gro+file_ndx+" "+folder_fill+"$fill_num\n"
            out_string += "cp "+folder_nvt+sim_nvt+".gro "+folder_gro+file_box+"\n"
            out_string += "cp "+folder_top+file_t_b+" "+folder_top+file_top+"\n"
            out_string += "mv "+folder_min+" "+folder_fill+"$fill_num\n"
            out_string += "mv "+folder_nvt+" "+folder_fill+"$fill_num\n"
            out_string += "mkdir "+folder_min+"\n"
            out_string += "mkdir "+folder_nvt+"\n"
            out_string += "cp "+folder_fill+"$fill_num/"+sim_min+"/"+sim_min+".job"+" "+folder_min+"\n"
            out_string += "cp "+folder_fill+"$fill_num/"+sim_nvt+"/"+sim_nvt+".job"+" "+folder_nvt+"\n"
            file_out.write(out_string)
            file_out.write("echo \"System "+self._box_link+" - Backed up equilibration ...\"\n\n")

            # Backup Analysis in case of automatic density
            if not all(x is None for x in [self._mols[mol][2] for mol in self._mols]):
                file_out.write("# Backup Analysis\n")
                out_string = "mv "+folder_ana+" "+folder_fill+"$fill_num\n"
                out_string += "mkdir "+folder_ana+"\n"
                out_string += "cp "+folder_fill+"$fill_num/ana/ana.* "+folder_ana+"\n"
                file_out.write(out_string)
                file_out.write("echo \"System "+self._box_link+" - Backed up analysis ...\"\n\n")

            # Fill box
            file_out.write("# Refill Box\n")
            for mol in self._mols:
                if self._mols[mol][0]=="fill":
                    out_string = "gmx_mpi insert-molecules "
                    out_string += "-f "+folder_gro+file_box+" "
                    out_string += "-o "+folder_gro+file_box+" "
                    out_string += "-ci "+folder_gro+self._struct[mol].split("/")[-1]+" "
                    out_string += "-try 10000 "
                    out_string += "-scale 0.57 "
                    out_string += "-nmol "
                    out_string += str(10000) if self._mols[mol][2] is None else "FILLDENS"
                    out_string += " >> logging.log 2>&1\n"
                    out_string += "python empty_grid.py "+folder_gro+" "+mol+" "+str(self._mols[mol][1])+"\n"
                    file_out.write(out_string)
            file_out.write("python sort.py "+folder_gro+"\n")
            file_out.write("echo \"System "+self._box_link+" - Refilled simulation box ...\"\n\n")

            # Update topology and create index
            self._topol_index(file_out, "../")
            file_out.write("\n")

            # Remove log and backup
            file_out.write("# Remove logs\n")
            file_out.write("rm "+folder_gro+"*#\n")
            file_out.write("rm logging.log\n\n")

            # Step a folder number forward
            file_out.write("# Step fill folder number\n")
            file_out.write("cp fill.sh temp.sh\n")
            file_out.write("sed -i \"s/fill_num=$fill_num/fill_num=$((fill_num+1))/\" temp.sh\n")
            if not all(x is None for x in [self._mols[mol][2] for mol in self._mols]):
                file_out.write("sed -i \"s/fill_num=$fill_num/fill_num=$((fill_num+1))/\" fillBackup.sh\n")
            file_out.write("mv temp.sh fill.sh\n")


    ##################
    # Public Methods #
    ##################
    def generate_files(self):
        """Generate structure files and shells.
        """
        # Create construction shell file
        self._structure()

        # Create structure folder
        utils.mkdirp(self._box_path+"_gro")

        # Copy structure files
        for mol in self._struct:
            file_link = self._struct[mol]
            if mol=="BOX":
                utils.copy(file_link, self._box_path+"_gro/"+"box.gro")
            elif mol=="GENERATE":
                utils.copy(file_link, self._box_path+"_gro/"+"generate.sh")
            elif mol=="PLUMED":
                utils.copy(file_link, self._box_path+"_gro/"+"plumed.dat")
            else:
                utils.copy(file_link, self._box_path+"_gro/"+file_link.split("/")[-1])

        # Pore simulation that needs to be filled
        if "fill" in [self._mols[mol][0] for mol in self._mols]:
            # Create filling backup folder
            utils.mkdirp(self._box_path+"_fill")

            # Create shell files
            self._fill()

            # Create fill backup for automatic filling
            if not all(x is None for x in [self._mols[mol][2] for mol in self._mols]):
                utils.copy(self._box_path+"_fill/fill.sh", self._box_path+"_fill/fillBackup.sh")

            # Copy empty script
            utils.copy(os.path.split(__file__)[0]+"/templates/empty_grid.py", self._box_path+"_fill/"+"empty_grid.py")
            utils.copy(os.path.split(__file__)[0]+"/templates/sort.py", self._box_path+"_fill/"+"sort.py")
