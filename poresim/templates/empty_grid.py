import sys
import yaml
import porems as pms


if __name__ == "__main__":
    # Get user input
    link = sys.argv[1]
    link = link+"/" if not link[-1] == "/" else link
    mol_name = sys.argv[2]
    mol_len = int(sys.argv[3])

    # Load pores object
    with open(link+"pore.yml", "r") as stream:
        pore = yaml.load(stream, Loader=yaml.FullLoader)

    # Get Diameter
    diam = pore["diameter"]
    focal = pore["centroid"]

    # Calculate pore position
    block = pms.utils.load(link+"pore.obj")
    mol_list = sum([block.get_mol_dict()[x] for x in sorted(list(block.get_mol_dict().keys()))], [])

    z_min = 1000000
    z_max = 0
    grid_mols = ["om", "ox", "si"]
    for mol in mol_list:
        if mol.get_name() in grid_mols:
            if mol.pos(0)[2] < z_min:
                z_min = mol.pos(0)[2]
            if mol.pos(0)[2] > z_max:
                z_max = mol.pos(0)[2]

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

    # Create molecule id lists
    mols = []
    counter = 0
    for atom in data:
        if mol_name in atom[1]:
            if counter == 0:
                mols.append([])
                counter = mol_len
            mols[-1].append(atom)
            counter -= 1

    # Find molecules stuck in grid
    del_list = []
    for i in range(len(mols)):
        for atom in mols[i]:
            # Calculate radial distance
            radial = pms.geom.length(pms.geom.vector([focal[0], focal[1], atom[6]], atom[4:]))

            # Check if position is within grid
            if atom[6] > z_min and atom[6] < z_max and radial > diam/2:
                del_list.append(i)
                break

    # Remove molecules in grid
    for i in sorted(del_list, reverse=True):
        mols.pop(i)

    # Feed new molecule list into data matrix
    for i in reversed(range(len(data))):
        if mol_name in data[i][1]:
            data.pop(i)
    for mol in mols:
        for atom in mol:
            data.append(atom)

    # Create new gro file
    with open(link+"box.gro", "w") as file_out:
        # Set title
        file_out.write(title)

        # Number of atoms
        file_out.write("%i" % len(data)+"\n")

        # Atoms
        num_a = 0
        num_m = 0
        last_m = -1
        for atom in data:
            # Set ids
            if not atom[0] == last_m:
                last_m = atom[0]
                num_m = num_m+1 if num_m < 99999 else 0
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
