################################################################################
# Actuate Class                                                                #
#                                                                              #
"""All necessary function for running the simulation."""
################################################################################


from re import A
import poresim.utils as utils


class Actuate:
    """This class creates needed shell files for running the simulation using
    GROMACS.

    Parameters
    ----------
    sim_link : string
        Simulation master folder link
    box_link : string
        Simulation box folder link
    cluster: dictionary
        Cluster information dictionary
    job : dictionary
        Job dictionary
    label : string
        Box label
    structure : dictionary
        Box structure dictionary
    """
    def __init__(self, sim_link, box_link, cluster, job, label, structure):
        # Initialize
        self._link = box_link
        self._sim_link = sim_link
        self._box_link = "./" if sim_link == box_link else "./"+box_link.split("/")[-2]+"/"
        self._clr_link = cluster["directory"]+sim_link.split("/")[-2]+"/"
        self._clr_link += "" if sim_link == box_link else box_link.split("/")[-2]+"/"
        self._cluster = cluster
        self._job = job
        self._label = label
        self._is_pore = "PORE" in structure
        self._is_plm = "PLUMED" in structure


    ###################
    # Private methods #
    ###################
    def _equilibration(self):
        """Create equilibration simulation shells for
        * **min** for energy minimization
        * **nvt** for temperature equilibration
        * **npt** for pressure equilibration
        """
        # Set simulation folder descriptors
        sim_min = "min"
        sim_nvt = "nvt"
        sim_npt = "npt"
        sim_all = {sim_min: sim_min, sim_nvt: sim_nvt, sim_npt: sim_npt}

        # Set folder names
        folder_gro = "../_gro/"
        folder_top = "../_top/"
        folder_mdp = "../_mdp/"

        # Set file names
        file_box = "box.gro"
        file_top = "topol.top"
        file_ndx = "index.ndx"

        # Get simulation properties
        np = {step: str(self._job[step]["np"]) for step in self._job}
        ntomp = {step: str(self._job[step]["ntomp"]) if "ntomp" in self._job[step] else "0" for step in self._job}
        nodes = {step: str(self._job[step]["nodes"]) for step in self._job}
        gpu = {step: self._job[step]["gpu"] if "gpu" in self._job[step] else False for step in self._job}
        wall = {step: self._job[step]["wall"] for step in self._job}

        # Set step dependent files
        conf = {"min": folder_gro+file_box,
                "nvt": "../"+sim_min+"/"+sim_min+".gro",
                "npt": "../"+sim_nvt+"/"+sim_nvt+".gro"}
        check = {"min": "",
                    "nvt": "",
                    "npt": "-t ../"+sim_nvt+"/"+sim_nvt+".cpt "}
        forward = {"min": sim_nvt, "nvt": sim_npt, "npt": ""}

        if "npt" not in self._job:
            forward["nvt"] = ""

        # Create job shells
        for step in self._job:
            if not step == "run":
                link_shell = self._link+step+"/"+step+".job"

                utils.mkdirp(self._link+step)
                utils.copy(self._job[step]["file"], link_shell)

                # Simulation options
                utils.replace(link_shell, "SIMULATIONNODES", nodes[step])
                utils.replace(link_shell, "SIMULATIONPROCS", np[step])
                utils.replace(link_shell, "SIMULATIONGPU", "gpu:" + str(gpu[step]) if gpu[step]!=False else "")
                utils.replace(link_shell, "SIMULATIONTIME", wall[step])
                utils.replace(link_shell, "SIMULATIONLABEL", self._label+"_"+step)
                utils.replace(link_shell, "COMMANDCHANGEDIR", "cd "+self._clr_link+step)

                # Grompp
                gro_string = ""
                gro_string += "gmx_mpi grompp "
                gro_string += "-f "+folder_mdp+sim_all[step]+".mdp "
                gro_string += "-c "+conf[step]+" "
                gro_string += check[step]
                gro_string += "-p "+folder_top+file_top+" "
                gro_string += "-n "+folder_gro+file_ndx+" " if self._is_pore else ""
                gro_string += "-po "+sim_all[step]+" "
                gro_string += "-o "+sim_all[step]
                gro_string += " -maxwarn 1" if self._is_pore else ""
                gro_string += "\n"

                # Mdrun
                gro_string += self._cluster["queuing"]["mpi"]+" "
                gro_string += str(int(np[step])*int(nodes[step]))+" " if self._cluster["queuing"]["add_np"] else ""
                gro_string += "gmx_mpi mdrun "
                gro_string += "-s "+sim_all[step]+".tpr "
                gro_string += "-ntomp "+ntomp[step]+" " if not ntomp[step]=="0" else ""
                gro_string += "-plumed ../_gro/plumed.dat " if self._is_plm else ""
                gro_string += "-deffnm "+sim_all[step]+" "
                gro_string += "-dlb no -c -e -g -o -cpo"

                # Insert into file
                utils.replace(link_shell, "COMMANDGROMACS", gro_string)

                # Add forward queuing
                with open(link_shell, "a") as file_out:
                    if not forward[step] == "":
                        file_out.write("\n")
                        file_out.write("sleep 20\n")
                        file_out.write("rm *#\n")
                        file_out.write("cd ../"+forward[step]+"/\n")
                        file_out.write(self._cluster["queuing"]["submit"]+" "+forward[step]+".job"+"\n")

        # Add to master shell
        with open(self._sim_link+"equilibrate.sh", "a") as file_out:
            file_out.write("cd "+self._box_link+sim_min+"/\n")
            file_out.write(self._cluster["queuing"]["submit"]+" "+sim_min+".job\n")
            back = "" if self._box_link == "./" else "../"
            file_out.write("cd "+back+"../\n")
            file_out.write("\n")

    def _simulation(self):
        """Create simulation job files. For the first run - for creating the tpr
        file - run0.job is used. For further simulations, job script run.job is
        created.

        Note that the latter automatically extends the current simulation.
        The breakpoint can be set using the *sim_num* variable in the shell
        file.
        """
        # Set simulation folder descriptor
        sim_nvt = "nvt"
        sim_npt = "npt"
        sim_run = "run"

        # Set folder names
        folder_gro = "../_gro/"
        folder_top = "../_top/"
        folder_mdp = "../_mdp/"

        # Set file names
        file_top = "topol.top"
        file_ndx = "index.ndx"

        # Get last step
        last = "../"+sim_npt+"/"+sim_npt if "npt" in self._job else "../"+sim_nvt+"/"+sim_nvt

        # Get simulation properties
        np = str(self._job["run"]["np"])
        nodes = str(self._job["run"]["nodes"])
        ntomp = str(self._job["run"]["ntomp"]) if "ntomp" in self._job["run"] else "0"
        gpu = self._job["run"]["gpu"] if "gpu" in self._job["run"] else False
        wall = self._job["run"]["wall"]

        # Create job shells
        for i in range(2):
            # Initialize
            is_first = i == 0
            shell_num = "0" if is_first else ""
            sim_num = i+1

            # Create shell
            link_shell = self._link+sim_run+"/"+sim_run+shell_num+".job"

            utils.mkdirp(self._link+sim_run)
            utils.copy(self._job["run"]["file"], link_shell)

            # Change variables
            utils.replace(link_shell, "SIMULATIONNODES", nodes)
            utils.replace(link_shell, "SIMULATIONPROCS", np)
            utils.replace(link_shell, "SIMULATIONGPU", ":gpus=1:exclusive_process" if gpu else "")
            utils.replace(link_shell, "SIMULATIONTIME", wall)
            utils.replace(link_shell, "SIMULATIONLABEL",self._label+"_"+str(sim_num))
            utils.replace(link_shell, "COMMANDCHANGEDIR", "cd "+self._clr_link+sim_run)

            # Grompp
            gro_string = ""
            if is_first:
                gro_string += "gmx_mpi grompp "
                gro_string += "-f "+folder_mdp+sim_run+".mdp "
                gro_string += "-c "+last+".gro "
                gro_string += "-t "+last+".cpt "
                gro_string += "-p "+folder_top+file_top+" "
                gro_string += "-n "+folder_gro+file_ndx+" " if self._is_pore else ""
                gro_string += "-po "+sim_run+" "
                gro_string += "-o "+sim_run+" "
                gro_string += " -maxwarn 1" if self._is_pore else ""
                gro_string += "\n"

            # Mdrun
            gro_string += self._cluster["queuing"]["mpi"]+" "
            gro_string += str(int(np)*int(nodes))+" " if self._cluster["queuing"]["add_np"] else ""
            gro_string += "gmx_mpi mdrun "
            gro_string += "-s "+sim_run+".tpr "
            gro_string += "-ntomp "+ntomp+" " if not ntomp=="0" else ""
            gro_string += "-plumed ../_gro/plumed.dat " if self._is_plm else ""
            gro_string += "-cpi "+sim_run+".cpt " if not is_first else ""
            gro_string += "-deffnm "+sim_run+" "
            gro_string += "-maxh "+str(self._job["run"]["maxh"])+" " if "maxh" in self._job["run"] else ""
            gro_string += "-dlb no -c -e -g -o -cpo"
            gro_string += " -append" if not is_first else ""

            # Insert into file
            utils.replace(link_shell, "COMMANDGROMACS", gro_string)

            # Add forward queuing
            with open(link_shell, "a") as file_out:
                file_out.write("\n")
                file_out.write("sim_num="+str(sim_num)+"\n")

                if not is_first:
                    file_out.write("\n")
                    file_out.write("cp "+sim_run+".job temp.job\n")
                    file_out.write("sed -i \"s/sim_num=$sim_num/sim_num=$((sim_num+1))/\" temp.job\n")
                    file_out.write("sed -i \"s/"+self._label+"_$sim_num/" +
                                  self._label+"_$((sim_num+1))/\" temp.job\n")
                    file_out.write("mv temp.job "+sim_run+".job\n")

                file_out.write("\n")
                file_out.write("sleep 20\n")
                out_string = "if "
                out_string += "[ -f "+sim_run+".cpt ] && "
                out_string += "[ $sim_num -lt "+str(self._job["run"]["runs"])+" ]; "
                out_string += "then "+self._cluster["queuing"]["submit"]+" "+sim_run+".job;"
                out_string += "fi\n"
                file_out.write(out_string)

        # Add to master shell
        with open(self._sim_link+"simulate.sh", "a") as file_out:
            file_out.write("cd "+self._box_link+sim_run+"/\n")
            file_out.write(self._cluster["queuing"]["submit"]+" "+sim_run+"0.job\n")
            back = "" if self._box_link == "./" else "../"
            file_out.write("cd "+back+"../\n")
            file_out.write("\n")

    def _analyze(self):
        """Create analyse job files for poreana.

        """
        # Get simulation properties
        np = str(int(20))
        nodes = str(int(1))
        wall = "10:00:00"

        # Create shell
        utils.mkdirp(self._link+"ana")
        link_shell = self._link+"ana/ana.job"
        utils.copy(self._job["run"]["file"], link_shell)

        # Change variables
        utils.replace(link_shell, "SIMULATIONNODES", nodes)
        utils.replace(link_shell, "SIMULATIONPROCS", np)
        utils.replace(link_shell, "SIMULATIONTIME", wall)
        utils.replace(link_shell, "SIMULATIONLABEL",self._label+"_"+str("ana"))
        utils.replace(link_shell, "COMMANDCHANGEDIR", "cd "+self._clr_link+str("ana"))

        # Insert into file
        utils.replace(link_shell, "COMMANDGROMACS", "python ana.py")

    ##################
    # Public Methods #
    ##################
    def generate_files(self):
        """Generate structure files and shells.
        """
        # Create equilibration files and folders
        self._equilibration()

        # Create simulation files for the production cycle
        if "run" in self._job:
            self._simulation()
        
            # Set ana.job file
            self._analyze()

