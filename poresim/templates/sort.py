import sys


if __name__ == "__main__":
    ### Anion ###
    # Set name of ions for ionic liquid simulations
    # an_ion = ""   # Set name "Im"
    # kat_ion = ""  # Set name "FAP"

    # Get user input
    link = sys.argv[1]
    link = link+"/" if not link[-1] == "/" else link

    # Read filled structure
    with open(link+"box.gro", "r") as file_in:
        counter = 0
        data = []
        for line in file_in:
            if counter < 2:
                if counter == 0:
                    title = line
                counter += 1
            elif not len(line.split()) == 3:
                atom = []
                atom.append(int(line[0:5]))
                atom.append(line[5:10])
                atom.append(line[10:15])
                atom.append(int(line[15:20]))
                atom.append(float(line[20:28]))
                atom.append(float(line[28:36]))
                atom.append(float(line[36:44]))

                data.append(atom)
            else:
                box = line

    # Create molecule lists
    mols = {}
    mol_names = []
    last_atom = None
    for atom in data:
        # Create new dictionary entry
        if not atom[1] in mols:
            mols[atom[1]] = []
            mol_names.append(atom[1])

        # Create new molecule
        if last_atom is None or not atom[0]==last_atom:
            mols[atom[1]].append([])

        # Add atom
        mols[atom[1]][-1].append(atom)
        last_atom = atom[0]

    ### Anion ###
    # Same amount of Kat-Ion and An-Ion for ionic liquid simulations
    # if an_ion is not None and kat_ion is not None:
    #     for mol_name in mol_names:
    #         if an_ion in mol_name:
    #             an_ion = mol_name
    #         if kat_ion in mol_name:
    #             kat_ion = mol_name
    #
    #     max_num_ions = min(len(mols[an_ion]), len(mols[kat_ion]))
    #
    #     mols[an_ion] = mols[an_ion][:max_num_ions]
    #     mols[kat_ion] = mols[kat_ion][:max_num_ions]

    # Create new gro file
    with open(link+"box.gro", "w") as file_out:
        # Set title
        file_out.write(title)

        # Number of atoms
        num_atoms = sum([sum([len(mol) for mol in mols[mol_name]]) for mol_name in mol_names])
        file_out.write("%i" % num_atoms+"\n")

        # Atoms
        num_a = 0
        num_m = 0
        for mol_name in mol_names:
            for mol in mols[mol_name]:
                num_m = num_m+1 if num_m < 99999 else 0
                for atom in mol:
                    # Set ids
                    num_a = num_a+1 if num_a < 99999 else 0

                    # Write atom string
                    out_string  = "%5i" % num_m  #  1- 5 (5)    Residue number
                    out_string += atom[1]        #  6-10 (5)    Residue short name
                    out_string += atom[2]        # 11-15 (5)    Atom name
                    out_string += "%5i" % num_a  # 16-20 (5)    Atom number
                    for j in range(3):           # 21-44 (3*8)  Coordinates
                        out_string += "%8.3f" % atom[4+j]

                    file_out.write(out_string+"\n")

        # Box
        file_out.write(box)
