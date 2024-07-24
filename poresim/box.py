################################################################################
# Box Class                                                                    #
#                                                                              #
"""Create simulation box object."""
################################################################################


class Box:
    """This class creates a hypothetical simulation box, to which molecules are
    added. Furthermore, all properties like temperature and specific job
    properties can be adjusted separate from the global simulation values.

    Parameters
    ----------
    name : string, optional
        Simulation name
    label : None, string, optional
        Job identification label
    """
    def __init__(self, name="box", label=None):
        # Initialize
        self._sim_dict = {}
        self._sim_dict["name"] = [name, name if label is None else label]

        self._sim_dict["struct"] = {}
        self._sim_dict["mols"] = {}
        self._sim_dict["topol"] = {"master":[], "top":[], "itp":[]}


    ##################
    # Public Methods #
    ##################
    def add_box(self, link):
        """Set simulation box structure link.

        Parameters
        ----------
        link : string
            Simulation box structure link
        """
        self._sim_dict["struct"]["BOX"] = link

    def add_pore(self, link):
        """Add pore object to structure folder.

        Parameters
        ----------
        link : string
            Pore object link
        """
        self._sim_dict["struct"]["PORE"] = link

    def add_mol(self, short, link, inp, num_atoms="gro", auto_dens=None, mass = None,  section="both",area= [], box = [], kwargs_gmx = {}):
        """Add a molecule to the simulation box. A unique short name and a
        structure-file link have to be given.

        Furthermore, the input is either a specific number of molecules or the
        input **fill** in order to continuously fill the simulation box with the
        specified molecule. This can be done automatically by defining the
        reference density to fill up to.

        The number atoms can either be specified numerically or extracted, which
        currently only works for the GROMACS file-type.

        Parameters
        ----------
        short : string
            Molecule short name as given in the topology
        link : string
            Molecule structure file link
        inp : integer, string
            Define number of molecules or **fill** for filling a pore or box
        num_atoms : integer, string, optional
            Number of atoms in the specified molecule or filetype for extraction
        auto_dens : float, None, optional
            Density in :math:`\\frac{\\text{kg}}{\\text{m}^3}`
        mass : float, None, optional
            molecular mass (:math:`\\frac{\\text{g}}{\\text{mol}}`) of the molecule to get an estimate of how many molecules you need to put in the box as the initial value for auto_dens.
        section : string, "both", optional
            'Both' the molecules in the pore and in the reservoir can be set. 
            If you want the molecules only in the pore, choose 'pore'. If you want them only in the reservoir, choose 'res'. 
            For a slit pore you can use wall to put the molecules as a layer on the SiO2 surface.
            Chose "box" for a box system.
        area : list, [], optional
            To set up a 2-phase box system without a pore, specify the areas of the section [a,b] where the molecules should be inserted using Gromacs.
        box : list, [], optional
            To set up a 2-phase system, specify the dimensions [x, y, z] of the simulation box.
        kwargs_gmx : dictonary, {}, optional
            dictonary to set up options for the gromacs "insert-molecules" command (adjust for example -try and -scale)
        """
        # Process input
        if not isinstance(inp, int) and not isinstance(inp, str):
            print("Wrong input for the molecule number ...")
            return

        if not isinstance(num_atoms, int) and not num_atoms in ["gro"]:
            print("Wrong input for the number of atoms ...")
            return

        if auto_dens is not None and not (isinstance(auto_dens, float) or isinstance(auto_dens, int)):
            print("Wrong input for the auto density option ...")
            return

        # Extract atom number
        if num_atoms == "gro":
            counter = 1
            with open(link, "r") as file_in:
                for line in file_in:
                    if counter == 2:
                        num_atoms = int(line)
                        break
                    else:
                        counter += 1

        # If Fill
        if inp=="fill" and not mass:
            print("If you use auto_dens you have to specifiy the molar mass of the molecule")
            return

        # Add to global list
        self._sim_dict["struct"][short] = link
        self._sim_dict["mols"][short] = [inp, num_atoms, auto_dens, mass,  section,area, box, kwargs_gmx]

    def add_struct(self, ident, link):
        """Add file link to the structure dictionary.

        Parameters
        ----------
        ident : string
            File identifier
        link : string
            File link
        """
        self._sim_dict["struct"][ident] = link

    def add_topol(self, link, top_type="itp"):
        """Add simulation topology files to simulation dictionary. For topology
        files, three different types exist with different handlings

        * **master** - Master topology file will be named topol.top with an additional backup
        * **top** - Topology file will be copied
        * **itp** - Topology file will be filtered for an itp file

        Parameters
        ----------
        link : string, list
            Link or a list of those containing a topology file link
        top_type : string
            Topology type
        """
        links = link if isinstance(link, list) else [link]
        for link in links:
            self._sim_dict["topol"][top_type].append(link)

    def add_charge_si(self, charge):
        """
        Add silicon charge for grid-atoms

        Parameters
        ----------
        charge : float
            Silicon charge for pore grid molecules
        """
        self._sim_dict["topol"]["charge"] = charge


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

    def set_name(self, name):
        """Set the simulation box name.

        Parameters
        ----------
        name : string
            Simulation box name
        """
        self._sim_dict["name"][0] = name

    def set_label(self, label):
        """Set the simulation box label, for job identification.

        Parameters
        ----------
        label : string
            Label name
        """
        self._sim_dict["name"][1] = label

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

    def get_name(self):
        """Return the simulation box name.

        Returns
        -------
        name : string
            Simulation box name
        """
        return self._sim_dict["name"][0]

    def get_label(self):
        """Return the simulation label.

        Returns
        -------
        label : string
            Simulation label
        """
        return self._sim_dict["name"][1]

    def get_job(self):
        """Return the dictionary containing all simulation information.

        Returns
        -------
        job : dictionary
            Simulation dictionary
        """
        return self._sim_dict["job"] if "job" in self._sim_dict.keys() else None

    def get_param(self):
        """Return the dictionary containing all systems parameters.

        Returns
        -------
        param : dictionary
            Parameter dictionary
        """
        return self._sim_dict["param"] if "param" in self._sim_dict.keys() else None

    def get_mols(self):
        """Return list of molecule numbers to be filled in the box.

        Returns
        -------
        mols : dictionary
            List of molecules
        """
        return self._sim_dict["mols"]

    def get_struct(self):
        """Return list of structure files.

        Returns
        -------
        struct : dictionary
            Molecule structure file links
        """
        return self._sim_dict["struct"]

    def get_topol(self):
        """Return list of topology files.

        Returns
        -------
        topol : list
            List of topology file links
        """
        return self._sim_dict["topol"]
