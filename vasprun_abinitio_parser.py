import numpy as np
import xmltodict

def ReadVasprun(filename):
	"""
	Input: path containing vasprun.xml file

	Output: Free energy for each abinitio MD runs, stress in GPa, Force in eV/Angstrom and Position of atoms

		Free energy: List containing free energy values in eV

		Stress: Output is a matrix where each row represent force for each atom in a crystal and column represent force in X, Y and Z direction

		Position: Output is a matrix where each row represent position of each atom and column represent X, Y and Z fractioal coordinates

		Lattice vector: Output is a matrix containing a1, a2 and a3 in each row

	"""
	structure_list=[]
	force_list=[]
	stress_list=[]
	free_energy_list=[]
	energy_list=[]
	lattice_vector_list=[]
	atom_info=[]    
	with open(filename) as file:
		# All of the abinitio MD data from VASP are in calculation section of vasprun.xml so that needs to be accessed first.
		doc=xmltodict.parse(file.read())
		atoms_info=[]        
		atom_des=doc['modeling']['atominfo']['array'][1]['set']['rc']
		atom_des_len=len(atom_des)
		for i in range(len(doc['modeling']['atominfo']['array'][1]['set']['rc'])):
			atoms_info.append(doc['modeling']['atominfo']['array'][1]['set']['rc'][i]['c'][0:2])            
		doc_calculation=doc['modeling']['calculation']        
		for i in range(len(doc_calculation)-1):
			tmp_structure=doc_calculation[i]['structure']['varray']['v']
			tmp_force=doc_calculation[i]['varray'][0]['v']
			tmp_stress=doc_calculation[i]['varray'][1]['v']
			tmp_lat_vec=doc_calculation[i]['structure']['crystal']['varray'][0]['v']
			free_energy_list.append(np.float64(doc_calculation[i]['energy']['i'][0]['#text'].strip()))
			energy_list.append(np.float64(doc_calculation[i]['energy']['i'][1]['#text'].strip()))
			structure_list.append(np.matrix(np.float64(list(map(lambda arg_string:arg_string.split(),tmp_structure)))))
			force_list.append(np.matrix(np.float64(list(map(lambda arg_string:arg_string.split(),tmp_force)))))
			stress_list.append(-0.0006241509125883258*np.matrix(np.float64(list(map(lambda arg_string:arg_string.split(),tmp_stress))))) # The number is conversion for stress in GPa, adapted from ASE
			lattice_vector_list.append(np.matrix(np.float64(list(map(lambda arg_string:arg_string.split(),tmp_lat_vec)))))         
	file.close()
	return lattice_vector_list, structure_list, force_list, stress_list, energy_list, free_energy_list, atoms_info
def cart_coordinate(output):
    result=[]
    cartesian=[]
    for i in range(len(output)):
        lat_vec=output[i][0]
        structure=output[i][1]
        new_cartesian=[]
        for j in range(len(structure)):
            cart_coord=[]
            for k in range(len(structure[j])):
                x=np.matmul(np.transpose(lat_vec[j][:,0]),np.transpose(structure[j][k,:]))[0][0,0]
                y=np.matmul(np.transpose(lat_vec[j][:,1]),np.transpose(structure[j][k,:]))[0][0,0]
                z=np.matmul(np.transpose(lat_vec[j][:,2]),np.transpose(structure[j][k,:]))[0][0,0]
                coords=[x,y,z]
                cart_coord.append(coords)
            new_cartesian.append(np.matrix(cart_coord))
        cartesian.append(new_cartesian)
    return cartesian
def writexyz(directory,arg_output,cartesian):
    os.mkdir(directory)
    dummy=0
    for i in range(len(arg_output)):
        for j in range(len(arg_output[i][0])):
            dummy=dummy+1
            with open(directory+"training_npt_new{}.xyz".format(dummy),"w") as file:
                file.write("{}".format(len(arg_output[i][1][j])))
                file.write("\n")
                a11=arg_output[i][0][j][0,0]
                a12=arg_output[i][0][j][0,1]
                a13=arg_output[i][0][j][0,2]
                a21=arg_output[i][0][j][1,0]
                a22=arg_output[i][0][j][1,1]
                a23=arg_output[i][0][j][1,2]
                a31=arg_output[i][0][j][2,0]
                a32=arg_output[i][0][j][2,1]
                a33=arg_output[i][0][j][2,2]
                energy=arg_output[i][4][j]
                free_energy=arg_output[i][5][j]
                sxx=output[i][3][j][0,0]
                syy=output[i][3][j][1,1]
                szz=output[i][3][j][2,2]
                syz=output[i][3][j][1,2]
                sxy=output[i][3][j][0,1]
                sxz=output[i][3][j][0,2]
                file.write('Lattice="{} {} {} {} {} {} {} {} {}" Properties=species:S:1:pos:R:3:forces:R:3 energy={} stress="{} {} {} {} {} {}" free_energy={} pbc="1 1 1"'.format(a11,a12,a13,a21,a22,a23,a31,a32,a33,energy,sxx,syy,szz,syz,sxy,sxz,free_energy))
                file.write("\n")
                for k in range(len(arg_output[i][1][j])):
                    x=cartesian[i][j][k,0]
                    y=cartesian[i][j][k,1]
                    z=cartesian[i][j][k,2]
                    fx=arg_output[i][2][j][k,0]
                    fy=arg_output[i][2][j][k,1]
                    fz=arg_output[i][2][j][k,2]
                    if k<int(arg_output[i][6][0][0]):
                        file.write("{} {} {} {} {} {} {}".format(arg_output[i][6][0][1],x,y,z,fx,fy,fz))
                    else:
                        file.write("{} {} {} {} {} {} {}".format(arg_output[i][6][1][1],x,y,z,fx,fy,fz))
                    file.write("\n")
            file.close()