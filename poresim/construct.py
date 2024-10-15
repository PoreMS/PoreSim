################################################################################
# Construct Class                                                              #
#                                                                              #
"""All necessary function for creating the finished simulation box."""
################################################################################


import os

import poresim.utils as utils
import numpy as np

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
        if "PORE" in struct:
            self._pore_props = utils.load(struct["PORE"])

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

    def _pos_dat(self):
        """
        Helper function to create position files for Gromacs to insert molecules in the reservoir, in the pore or in a specify area of a box system
        """
        for mol in self._mols:
            # Position file for reservoir
            if (self._mols[mol][0]=="fill" and not self._mols[mol][6] and self._mols[mol][4] in ["res", "both"]) and "PORE" in self._struct:
                num = int(self._mols[mol][2]/self._mols[mol][3]/10*6.022*self._pore_props["system"]["dimensions"][0]*self._pore_props["system"]["dimensions"][1]*self._pore_props["system"]["reservoir"]*2*0.5)
                with open(self._box_path +"_gro/" + "position_{}.dat".format(mol), "w") as file_out:
                    for i in range(int(num/2)):
                        out_string = str(self._pore_props["system"]["dimensions"][0]/2) + " "
                        out_string += str(self._pore_props["system"]["dimensions"][1]/2) + " "
                        out_string += str(self._pore_props["system"]["reservoir"]/2) + "\n"
                        file_out.write(out_string)
                    for i in range(int(num/2)):
                        out_string = str(self._pore_props["system"]["dimensions"][0]/2) + " "
                        out_string += str(self._pore_props["system"]["dimensions"][1]/2) + " "
                        out_string += str(self._pore_props["system"]["dimensions"][2] - self._pore_props["system"]["reservoir"]/2) + "\n"
                        file_out.write(out_string)
                    file_out.close()
                # Position files for pore
                for pore_id in self._pore_props.keys():
                    if pore_id[:5]=="shape":
                        if (self._pore_props[pore_id]["parameter"]["central"]==[0,0,1]) and (self._mols[mol][4] in ["pore", "both"]): 
                            num_pore = int(self._mols[mol][2]/self._mols[mol][3]/10*6.022*np.pi*self._pore_props[pore_id]["diameter"]**2/4*(self._pore_props[pore_id]["parameter"]["length"]))
                            with open(self._box_path +"_gro/" + "position_{}_{}.dat".format(pore_id,mol), "w") as file_out:
                                for i in range(num_pore):
                                    out_string = str(self._pore_props[pore_id]["parameter"]["centroid"][0]) + " "
                                    out_string += str(self._pore_props[pore_id]["parameter"]["centroid"][1]) + " "
                                    out_string += str(self._pore_props["system"]["reservoir"] + self._pore_props[pore_id]["parameter"]["centroid"][2]) + "\n"
                                    file_out.write(out_string)
                                file_out.close()
            elif (self._mols[mol][0]=="fill" and not self._mols[mol][6] and self._mols[mol][4] in ["pore"]) and "PORE" in self._struct:
                for pore_id in self._pore_props.keys():
                        if pore_id[:5]=="shape":
                            if (self._pore_props[pore_id]["parameter"]["central"]==[0,0,1]) and (self._mols[mol][4] in ["pore", "both"]): 
                                num_pore = int(self._mols[mol][2]/self._mols[mol][3]/10*6.022*np.pi*self._pore_props[pore_id]["diameter"]**2/4*(self._pore_props[pore_id]["parameter"]["length"]))
                                with open(self._box_path +"_gro/" + "position_{}_{}.dat".format(pore_id,mol), "w") as file_out:
                                    for i in range(num_pore):
                                        out_string = str(self._pore_props[pore_id]["parameter"]["centroid"][0]) + " "
                                        out_string += str(self._pore_props[pore_id]["parameter"]["centroid"][1]) + " "
                                        out_string += str(self._pore_props["system"]["reservoir"] + self._pore_props[pore_id]["parameter"]["centroid"][2]) + "\n"
                                        file_out.write(out_string)
                                    file_out.close()
            # Position file for put a specific number of molecules in pore or reservoir 
            elif self._mols[mol][0]!="fill" and not self._mols[mol][5] and not self._mols[mol][6] and "PORE" in self._struct:   
                j=0
                # Check how much molecules have to set in pore or in reservoir if a number is given
                for pore_id in self._pore_props.keys():
                        if pore_id[:5]=="shape":
                            j = j + 1 
                if (self._mols[mol][4]=="both"):
                    num = [int(self._mols[mol][0]/(j+2)) for i in range(j+2)]

                    if sum(num) != self._mols[mol][0]:
                        num[-1] = num[-1] + abs(self._mols[mol][0]-sum(num))

                elif (self._mols[mol][4]=="res"):
                    num = [int(self._mols[mol][0]/(2)) for i in range(2)]
                    if sum(num) != self._mols[mol][0]:
                        num[-1] = num[-1] + abs(self._mols[mol][0]-sum(num))

                elif (self._mols[mol][4]=="pore"):
                    num = [int(self._mols[mol][0]/(j)) for i in range(j)]
                    if sum(num) != self._mols[mol][0]:
                        num[-1] = num[-1] + abs(self._mols[mol][0]-sum(num))
                elif (self._mols[mol][4]=="box"):
                    num = [int(self._mols[mol][0]/(j)) for i in range(j)]
                    if sum(num) != self._mols[mol][0]:
                        num[-1] = num[-1] + abs(self._mols[mol][0]-sum(num))

                # Generate the position files for the number of molecules to set in the pores
                if self._mols[mol][4] in ["pore", "both"]:
                    for pore_id,j in zip(self._pore_props.keys(),range(j)):
                        if pore_id[:5]=="shape":
                            if self._pore_props[pore_id]["parameter"]["central"]==[0,0,1]:
                                with open(self._box_path +"_gro/" + "position_{}_{}.dat".format(pore_id,mol), "w") as file_out:
                                    for i in range(num[j]):
                                        out_string = str(self._pore_props[pore_id]["parameter"]["centroid"][0]) + " "
                                        out_string += str(self._pore_props[pore_id]["parameter"]["centroid"][1]) + " "
                                        out_string += str(self._pore_props["system"]["reservoir"] + self._pore_props[pore_id]["parameter"]["centroid"][2]) + "\n"
                                        file_out.write(out_string)
                                    file_out.close()

                # Generate the position files for the number of molecules to set in reservoir
                if self._mols[mol][4] in ["res", "both"]:
                    with open(self._box_path +"_gro/" + "position_{}.dat".format(mol), "w") as file_out:
                        for i in range(num[0]):
                            out_string = str(self._pore_props["system"]["dimensions"][0]/2) + " "
                            out_string += str(self._pore_props["system"]["dimensions"][1]/2) + " "
                            out_string += str(self._pore_props["system"]["reservoir"]/2) + "\n"
                            file_out.write(out_string)
                        for i in range(num[-1]):
                            out_string = str(self._pore_props["system"]["dimensions"][0]/2) + " "
                            out_string += str(self._pore_props["system"]["dimensions"][1]/2) + " "
                            out_string += str(self._pore_props["system"]["dimensions"][2] - self._pore_props["system"]["reservoir"]/2) + "\n"
                            file_out.write(out_string)
                        file_out.close()

                # If slit pore you can set the molecules on wall 
                elif self._mols[mol][4]=="wall":
                    for pore_id,j in zip(self._pore_props.keys(),range(j)):
                        if pore_id[:5]=="shape":
                                with open(self._box_path +"_gro/" + "position_{}_{}.dat".format(pore_id,mol), "w") as file_out:
                                    for i in range(int(self._mols[mol][0]/2)):
                                        out_string = str(self._pore_props[pore_id]["parameter"]["centroid"][0]) + " "
                                        out_string += str(self._pore_props[pore_id]["parameter"]["centroid"][1]-self._pore_props[pore_id]["parameter"]["height"]/2*0.75) + " "
                                        out_string += str(self._pore_props["system"]["reservoir"] + self._pore_props[pore_id]["parameter"]["centroid"][2]) + "\n"
                                        file_out.write(out_string)
                                    for i in range(int(self._mols[mol][0]/2)):
                                        out_string = str(self._pore_props[pore_id]["parameter"]["centroid"][0]) + " "
                                        out_string += str(self._pore_props[pore_id]["parameter"]["centroid"][1]+self._pore_props[pore_id]["parameter"]["height"]/2*0.75) + " "
                                        out_string += str(self._pore_props["system"]["reservoir"] + self._pore_props[pore_id]["parameter"]["centroid"][2]) + "\n"
                                        file_out.write(out_string)

            # Position file for box system 
            elif self._mols[mol][0]=="fill" and self._mols[mol][5] and not self._mols[mol][4]=="wall":
                for area,i in zip(self._mols[mol][5], range(len(self._mols[mol][5]))):
                    num = int((self._mols[mol][2]/self._mols[mol][3]/10*6.022*self._mols[mol][6][0]*self._mols[mol][6][1]*(area[1]-area[0]))*0.8)
                    with open(self._box_path +"_gro/" + "position_{}_area{}.dat".format(mol,i), "w") as file_out:
                        for i in range(num):
                            out_string = str(self._mols[mol][6][0]/2) + " "
                            out_string += str(self._mols[mol][6][1]/2) + " "
                            out_string += str((area[1]+area[0])/2) + "\n"
                            file_out.write(out_string)
                        file_out.close()
            

            elif self._mols[mol][0]=="fill" and not self._mols[mol][5] and not self._mols[mol][4]=="wall":
                num = int((self._mols[mol][2]/self._mols[mol][3]/10*6.022*self._mols[mol][6][0]*self._mols[mol][6][1]*self._mols[mol][6][2])*0.8)
                with open(self._box_path +"_gro/" + "position_{}.dat".format(mol), "w") as file_out:
                    for i in range(num):
                        out_string = str(self._mols[mol][6][0]/2) + " "
                        out_string += str(self._mols[mol][6][1]/2) + " "
                        out_string += str((self._mols[mol][6][2])/2) + "\n"
                        file_out.write(out_string)
                file_out.close()
            
            elif self._mols[mol][0]!="fill" and not self._mols[mol][5]:
                if not self._mols[mol][6]:
                    print("If you fill a system with on molecule specify in add_mol for every molecule the box dimension")
                    exit()
                num = self._mols[mol][0]
                with open(self._box_path +"_gro/" + "position_{}.dat".format(mol), "w") as file_out:
                    for i in range(num):
                        out_string = str(self._mols[mol][6][0]/2) + " "
                        out_string += str(self._mols[mol][6][1]/2) + " "
                        out_string += str((self._mols[mol][6][2])/2) + "\n"
                        file_out.write(out_string)
                file_out.close()

            # Position file for box system to put molecules in a certain area
            elif self._mols[mol][0]!="fill" and self._mols[mol][5]:
                num = self._mols[mol][0]
                for area,i in zip(self._mols[mol][5], range(len(self._mols[mol][5]))):
                    with open(self._box_path +"_gro/" + "position_{}_area{}.dat".format(mol,i), "w") as file_out:
                        for i in range(int(num/len(self._mols[mol][5]))):
                            out_string = str(self._mols[mol][6][0]/2) + " "
                            out_string += str(self._mols[mol][6][1]/2) + " "
                            out_string += str((area[1]+area[0])/2) + "\n"
                            file_out.write(out_string)
                    file_out.close()


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
            
            # First lines for standard gromacs command
            gmx_standard =  "gmx_mpi insert-molecules " + "-f "+folder_gro+file_box+" " +  "-o "+folder_gro+file_box+" "
            
            # Loop over the molecules which you want to insert into the system
            for mol in self._mols:
                file_out.write("\n####### " + mol + " #########\n")
                # Fill box
                if "PORE" in self._struct:
                    if (self._mols[mol][0]=="fill") and (self._mols[mol][3]!=None):
                        num = int(self._mols[mol][2]/self._mols[mol][3]/10*6.022*self._pore_props["system"]["dimensions"][0]*self._pore_props["system"]["dimensions"][0]*self._pore_props["system"]["reservoir"]*2)
                    else:
                        num = 0
                    if (self._mols[mol][4] in ["res","both"]):
                        # Fill reservoir 
                        file_out.write("# Fill Reservoir " + mol +" \n")
                        out_string = gmx_standard
                        out_string += "-ci "+folder_gro+self._struct[mol].split("/")[-1]+" "
                        out_string += "-dr "+ str(self._pore_props["system"]["dimensions"][0]/2) + " " + str(self._pore_props["system"]["dimensions"][1]/2) + " " + str(self._pore_props["system"]["reservoir"]/2) + " "
                        out_string += "-ip "+ folder_gro + "position_{}.dat".format(mol)  +" "
                        out_string += "-nmol "+str(int(self._mols[mol][0])) if not self._mols[mol][0]=="fill" else "-nmol "+str(num) 
                        for key,value in self._mols[mol][-1].items():
                            out_string +=  " "
                            out_string += key + " " + str(value) + " "
                        out_string += " >> logging.log 2>&1\n"
                        file_out.write(out_string)
                        file_out.write("echo \"Filled reservoir with " + mol + " ...\"\n\n")

                    # Fill pore area
                    if self._mols[mol][4]=="pore" or self._mols[mol][4]=="both":
                        file_out.write("# Fill Pore " + mol +" \n")
                        for pore_id in self._pore_props.keys():
                            if pore_id[:5]=="shape":
                                if self._pore_props[pore_id]["parameter"]["central"]==[0,0,1]:
                                    if (self._mols[mol][0]=="fill" and self._mols[mol][3]!=None ):
                                        if self._pore_props[pore_id]["shape"]=="SLIT":
                                            num_pore = int(self._mols[mol][2]/self._mols[mol][3]/10*6.022*self._pore_props["system"]["dimensions"][0]*self._pore_props[pore_id]["diameter"]*(self._pore_props[pore_id]["parameter"]["length"]))
                                        else:
                                            num_pore = int(self._mols[mol][2]/self._mols[mol][3]/10*6.022*np.pi*self._pore_props[pore_id]["diameter"]**2/4*(self._pore_props[pore_id]["parameter"]["length"]))
                                    else:
                                        num_pore = 0
                                    if self._mols[mol][4]!="res":
                                        out_string = gmx_standard
                                        out_string += "-ci "+folder_gro+self._struct[mol].split("/")[-1]+" "
                                        if self._pore_props[pore_id]["shape"]=="SLIT":
                                            out_string += "-dr "+ str(0.90*self._pore_props["system"]["dimensions"][0]/2) + " " + str(0.90*self._pore_props[pore_id]["diameter"]/2) + " " + str((0.9*self._pore_props[pore_id]["parameter"]["length"])/2) + " "
                                        else:
                                            out_string += "-dr "+ str(0.50*np.sqrt(0.9 *self._pore_props[pore_id]["diameter"]**2)/2) + " " + str(0.50*np.sqrt(0.9 *self._pore_props[pore_id]["diameter"]**2)/2) + " " + str((0.9*self._pore_props[pore_id]["parameter"]["length"])/2) + " "
                                        out_string += "-ip "+ folder_gro + "position_{}_{}.dat".format(pore_id,mol)  +" " 
                                        out_string += "-nmol "+str(int(self._mols[mol][0])) if not self._mols[mol][0]=="fill" else "-nmol "+str(num_pore) + "  "
                                        for key,value in self._mols[mol][-1].items():
                                            out_string +=  " "
                                            out_string += key + "  " + str(value) + " "
                                        out_string += " >> logging.log 2>&1\n"
                                        file_out.write(out_string)
                        file_out.write("echo \"Filled pore with " + mol + " ...\"\n\n")
                    
                    # Fill wall
                    elif self._mols[mol][4] == "wall":
                        file_out.write("# Fill Pore Wall " + mol +" \n")
                        for pore_id in self._pore_props.keys():
                            if pore_id[:5]=="shape":
                                if self._pore_props[pore_id]["parameter"]["central"]==[0,0,1]:
                                    out_string = gmx_standard
                                    out_string += "-ci "+folder_gro+self._struct[mol].split("/")[-1]+" "
                                    out_string += "-dr "+ str(0.9*self._pore_props["system"]["dimensions"][0]/2) + " " + str(self._pore_props[pore_id]["diameter"]/2*0.1) + " " + str((self._pore_props[pore_id]["parameter"]["length"])/2) + " "
                                    out_string += "-ip "+ folder_gro + "position_{}_{}.dat".format(pore_id,mol)  +" " 
                                    out_string += "-nmol "+str(int(self._mols[mol][0])) if not self._mols[mol][0]=="fill" else "-nmol "+str(num_pore) + "  "
                                    for key,value in self._mols[mol][-1].items():
                                        out_string +=  " "
                                        out_string += key + "  " + str(value) + " "
                                    out_string += " >> logging.log 2>&1\n"
                                    file_out.write(out_string)
                                file_out.write("echo \"Filled " + pore_id + " " +  mol + " ...\"\n\n")
                # If only box system
                else:
                    # Fill box system in a certain area
                    file_out.write("# Fill Box\n")
                    if self._mols[mol][5]:
                        for area,i in zip(self._mols[mol][5], range(len(self._mols[mol][5]))):
                            if self._mols[mol][0]=="fill":
                                num = int(self._mols[mol][2]/self._mols[mol][3]/10*6.022*self._mols[mol][6][0]*self._mols[mol][6][1]*(area[1]-area[0]))
                            else:
                                num = self._mols[mol][0]
                            out_string = gmx_standard
                            out_string += "-ci "+folder_gro+self._struct[mol].split("/")[-1]+" "
                            out_string += "-dr "+ str(self._mols[mol][6][0]/2) + " " + str(self._mols[mol][6][1]/2) + " " + str((area[1]-area[0])/2) + " "
                            out_string += "-ip "+ folder_gro + "position_{}_area{}.dat".format(mol,i)  +" " 
                            out_string += "-nmol "+str(int(self._mols[mol][0])) if not self._mols[mol][0]=="fill" else "-nmol "+str(num) + " "
                            for key,value in self._mols[mol][-1].items():
                                out_string += " "
                                out_string += key + " " + str(value) + " "
                            out_string += " >> logging.log 2>&1\n"
                            file_out.write(out_string)   
                    else: 
                        # Fill box with density or specific number of molecules
                        if self._mols[mol][0]=="fill":
                            num = int(self._mols[mol][2]/self._mols[mol][3]/10*6.022*self._mols[mol][6][0]*self._mols[mol][6][1]*self._mols[mol][6][2])
                        else:
                            num = self._mols[mol][0]
                        out_string = gmx_standard
                        out_string += "-ci "+folder_gro+self._struct[mol].split("/")[-1]+" "
                        if self._mols[mol][0]=="fill":
                            out_string += "-dr "+ str(self._mols[mol][6][0]/2) + " " + str(self._mols[mol][6][1]/2) + " " + str((self._mols[mol][6][2])/2) + " "
                            out_string += "-ip "+ folder_gro + "position_{}.dat".format(mol)  +" " 
                        out_string += "-nmol "+str(int(num))
                        for key,value in self._mols[mol][-1].items():
                            out_string += " "
                            out_string += key + " " + str(value) + " "
                        out_string += " >> logging.log 2>&1\n"
                        file_out.write(out_string)

            #if "fill" in [self._mols[mol][0] for mol in self._mols]:
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
            if "PORE" in self._struct:
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
                if self._mols[mol][0]=="fill" and self._mols[mol][5]:
                    for area,i in zip(self._mols[mol][5], range(len(self._mols[mol][5]))):
                        out_string = "gmx_mpi insert-molecules "
                        out_string += "-f "+folder_gro+file_box+" "
                        out_string += "-o "+folder_gro+file_box+" "
                        out_string += "-ci "+folder_gro+self._struct[mol].split("/")[-1]+" "
                        out_string += "-try 1000 "
                        out_string += "-scale 0.47 "
                        out_string += "-dr "+ str(self._mols[mol][6][0]/2) + " " + str(self._mols[mol][6][1]/2) + " " + str((area[1]-area[0])/2) + " "
                        out_string += "-ip "+ folder_gro + "position_{}_area{}.dat".format(mol,i)  +" " 
                        out_string += "-nmol "
                        out_string += str(10000) if self._mols[mol][2] is None else ("FILLDENS_" + mol)
                        out_string += " >> logging.log 2>&1\n" 
                        file_out.write(out_string)
                elif self._mols[mol][0]=="fill" and not self._mols[mol][5]:
                    out_string = "gmx_mpi insert-molecules "
                    out_string += "-f "+folder_gro+file_box+" "
                    out_string += "-o "+folder_gro+file_box+" "
                    out_string += "-ci "+folder_gro+self._struct[mol].split("/")[-1]+" "
                    out_string += "-try 1000 "
                    out_string += "-scale 0.47 "
                    if "PORE" in self._struct:
                        out_string += "-dr "+ str(self._pore_props["system"]["dimensions"][0]/2) + " " + str(self._pore_props["system"]["dimensions"][1]/2) + " " + str(self._pore_props["system"]["reservoir"]/2) + " "
                    else:
                        out_string += "-dr "+ str(self._mols[mol][6][0]/2) + " " + str(self._mols[mol][6][1]/2) + " " + str(self._mols[mol][6][2]/2) + " "
                    out_string += "-ip "+ folder_gro + "position_{}.dat".format(mol)  +" " 
                    out_string += "-nmol "
                    out_string += 0 if self._mols[mol][2] is None else ("FILLDENS_" + mol)
                    out_string += " >> logging.log 2>&1\n"
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
  
        # If molecules have to put in a specific area of a pore system 
        if "PORE" in self._struct:
            if ("wall" or "pore" in [self._mols[mol][4] for mol in self._mols]):
                if not "box" in [self._mols[mol][4] for mol in self._mols]:
                    self._pos_dat()

        # If molecules have to set in a specific area of a box system
        for mol in self._mols:
            if self._mols[mol][5]:
                self._pos_dat()

        # Pore simulation that needs to be filled
        if "fill" in [self._mols[mol][0] for mol in self._mols]:
            # Create filling backup folder
            utils.mkdirp(self._box_path+"_fill")

            # Create position file
            self._pos_dat()

            # Create shell files
            self._fill()

            # Create fill backup for automatic filling
            if not all(x is None for x in [self._mols[mol][2] for mol in self._mols]):
                utils.copy(self._box_path+"_fill/fill.sh", self._box_path+"_fill/fillBackup.sh")

            utils.copy(os.path.split(__file__)[0]+"/templates/sort.py", self._box_path+"_fill/"+"sort.py")
