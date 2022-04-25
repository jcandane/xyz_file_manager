import numpy as np 

π = np.pi 

def Bohr(x): ### given Ångström get Bohr
    return 1.889726125 * x 

def Ångström(x): ### given Bohr get Ångström
    return x / 1.889726125

def print2darray_tostring(Array_2D, ZZ):
    """ ZZ is the input from UHF object (LiH.Z) """
    output = ''
    for index, atom in enumerate(Array_2D):
        output += "  " + ZZ[index]
        for xyz in atom:
            number = np.format_float_positional(np.float32(xyz), unique=False, precision=6, pad_left = True)
            output += "\t  " + number
        output += ";\n"
    return output.rstrip()

def xyz_reader(filename, input_unit="Å", output_unit="å", mass=False, names=False): #currentlogfile = 'Ar_1000.xyz'#'Ne_1000.xyz'#'Na_7.xyz'#'Na_1000.xyz'# #currentlogfile = 'Xe_1Million.xyz' #currentlogfile = 'ArTest20.xyz' #ArTest.xyz
    """
    GIVEN:  .xyz file
            **input  unit
            **output unit
            **mass (Boolean, in a.u. units)
            **name (Boolean, element name)
    GET:    Z (atomic number)
            R_ix (position np.array, in chosen units)
            m_i (masses)
            names (atomic names on periodic table)
    """

    units = {"Å":0, "Ångström":0, "Angstrom":0, "angstrom":0, "ANGSTROM":0, "ÅNGSTRÖM":0, 
             "nm":1, "NM":1, "nanometer":1, "NANOMETER":1, 
             "å":2, "bohr":2, "Bohr":2, "BOHR":2, "a0":2}

    unitconversionmatrix = np.array([[1., 0.1, 1.8897259885], [10., 1., 18.8972598858], [1/1.8897259885,1/18.8972598858,1.]])
    u = unitconversionmatrix[ units[input_unit.strip()], units[output_unit.strip()] ]

    Element_Names = {'e ':0,  0:'e ',  'H ':  1,  1:'H ',  'He':  2,  2:'He',
                  'Li':  3,  3:'Li',  'Be':  4,  4:'Be',  'B ':  5,  5:'B ',  'C ':  6,  6:'C ',  'N ':  7,  7:'N ',  'O ':  8,  8:'O ',  'F ':  9,  9:'F ',  'Ne': 10, 10:'Ne', 
                  'Na': 11, 11:'Na',  'Mg': 12, 12:'Mg',  'Al': 13, 13:'Al',  'Si': 14, 14:'Si',  'P ': 15, 15:'P ',  'S ': 16, 16:'S ',  'Cl': 17, 17:'Cl',  'Ar': 18, 18:'Ar',
                  'K ': 19, 19:'K ',  'Ca': 20, 20:'Ca',  'Sc': 21, 21:'Sc',  'Ti': 22, 22:'Ti',  'V ': 23, 23:'V ',  'Cr': 24, 24:'Cr',  'Mn': 25, 25:'Mn',  'Fe': 26, 26:'Fe',  'Co': 27, 27:'Co',  'Ni': 28, 28:'Ni',  'Cu': 29, 29:'Cu',  'Zn': 30, 30:'Zn',  'Ga': 31, 31:'Ga',  'Ge': 32, 32:'Ge',  'As': 33, 33:'As',  'Se': 34, 34:'Se',  'Br': 35, 35:'Br',  'Kr': 36, 36:'Kr', 
                  'Rb': 37, 37:'Rb',  'Sr': 38, 38:'Sr',  'Y ': 39, 39:'Y ',  'Zr': 40, 40:'Zr',  'Nb': 41, 41:'Nb',  'Mo': 42, 42:'Mo',  'Tc': 43, 43:'Tc',  'Ru': 44, 44:'Ru',  'Rh': 45, 45:'Rh',  'Pd': 46, 46:'Pd',  'Ag': 47, 47:'Ag',  'Cd': 48, 48:'Cd',  'In': 49, 49:'In',  'Sn': 50, 50:'Sn',  'Sb': 51, 51:'Sb',  'Te': 52, 52:'Te',  'I ': 53, 53:'I ',  'Xe': 54, 54:'Xe',
                  'Cs': 55, 55:'Cs',  'Ba': 56, 56:'Ba',  'La': 57, 57:'La',  'Ce': 58, 58:'Ce',  'Pr': 59, 59:'Pr',  'Nd': 60, 60:'Nd',  'Pm': 61, 61:'Pm',  'Sm': 62, 62:'Sm',  'Eu': 63, 63:'Eu',  'Gd': 64, 64:'Gd',  'Tb': 65, 65:'Tb',  'Dy': 66, 66:'Dy',  'Ho': 67, 67:'Ho',  'Er': 68, 68:'Er',  'Tm': 69, 69:'Tm',  'Yb': 70, 70:'Yb',  'Lu': 71, 71:'Lu',  'Hf': 72, 72:'Hf',  'Ta': 73, 73:'Ta',  'W ': 74, 74:'W ',  'Re': 75, 75:'Re',  'Os': 76, 76:'Os',  'Ir': 77, 77:'Ir',  'Pt': 78, 78:'Pt',  'Au': 79, 79:'Au', 'Hg': 80, 80:'Hg',  'Tl': 81, 81:'Tl',  'Pb': 82, 82:'Pb',  'Bi': 83, 83:'Bi',  'Po': 84, 84:'Po', 'At': 85, 85:'At', 'Rn': 86, 86:'Rn', 
                  'Fr': 87, 87:'Fr',  'Ra': 88, 88:'Ra',  'Ac': 89, 89:'Ac',  'Th': 90, 90:'Th',  'Pa': 91, 91:'Pa',  'U ': 92, 92:'U ',  'Np': 93, 93:'Np',  'Pu': 94, 94:'Pu',  'Am': 95, 95:'Am',  'Cm': 96, 96:'Cm',  'Bk': 97, 97:'Bk',  'Cf': 98, 98:'Cf',  'Es': 99, 99:'Es',  'Fm':100,100:'Fm',  'Md':101,101:'Md',  'No':102,102:'No',  'Lr':103,103:'Lr',  'Rf':104,104:'Rf',  'Db':105,105:'Db',  'Sg':106,106:'Sg',  'Bh':107,107:'Bh',  'Hs':108,108:'Hs',  'Mt':109,109:'Mt',  'Ds':110,110:'Ds',  'Rg':111,111:'Rg', 'Cn':112,112:'Cn',  'Nh':113,113:'Nh',  'Fl':114,114:'Fl',  'Mc':115,115:'Mc',  'Lv':116,116:'Lv', 'Ts':117,117:'Ts', 'Og':118,118:'Og'}

    Atom_Names = np.array(["e ", "H ", "He", 
    "Li", "Be", "B ", "C ", "N ", "O ", "F ", "Ne", 
    "Na", "Mg", "Al", "Si", "P ", "S ", "Cl", "Ar", 
    "K ", "Ca", "Sc", "Ti", "V ", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", 
    "Rb", "Sr", "Y ", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I ", "Xe", 
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W ", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", 
    "Fr", "Ra", "Ac", "Th", "Pa", "U ", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"])

    Z_mass = np.array([ 1.,   1837.,   7297.,  
                   12650., 16427.,  19705.,  21894.,  25533.,  29164.,  34631., 36785.,  
                   41908., 44305.,  49185.,  51195.,  56462.,  58441.,  64621.,  72820.,  
                   71271., 73057.,  81949.,  87256.,  92861.,  94782., 100145., 101799., 107428., 106990., 115837., 119180., 127097., 132396., 136574., 143955., 145656., 152754.,
                   155798., 159721., 162065., 166291., 169357., 174906., 176820., 184239., 187586., 193991., 196631., 204918., 209300., 216395., 221954., 232600., 231331., 239332., 
                   242270., 250331., 253208., 255415., 256859., 262937., 264318., 274089., 277013., 286649., 289702., 296219., 300649., 304894., 307947., 315441., 318945., 325367., 329848., 335119., 339434., 346768., 350390., 355616., 359048., 365656., 372561., 377702., 380947., 380983., 382806., 404681.,
                   406504., 411972., 413795., 422979., 421152., 433900., 432024., 444784., 442961., 450253., 450253., 457545., 459367., 468482., 470305., 472128., 477596., 486711., 492179., 490357., 492179., 492179., 506763., 512231., 512231., 519523., 521346., 526814., 526814., 534106., 534106., 535929.])

    file = open(filename,'r')
    lines = file.readlines()
    lines.pop(0)
    lines.pop(0)
    
    Z = np.array([], dtype=np.int8)
    x = np.array([])
    y = np.array([])
    z = np.array([])

    for element in lines:
        a_line_in_lines = element.split()
    
        Z = np.append(Z , [np.int8(Element_Names[element[0] + element[1]])], axis = 0) 
        x = np.append(x, [(float(a_line_in_lines[1]))], axis = 0)
        y = np.append(y, [(float(a_line_in_lines[2]))], axis = 0)
        z = np.append(z, [(float(a_line_in_lines[3]))], axis = 0)
    file.close()

    if mass == True:
        if names == True:
            return Z, u * np.stack((x, y, z)).T, Z_mass[Z], Atom_Names[Z]
        else:
            return Z, u * np.stack((x, y, z)).T, Z_mass[Z]
    else:
        if names == True:
            return Z, u * np.stack((x, y, z)).T, Atom_Names[Z]
        else:
            return Z, u * np.stack((x, y, z)).T

def Get_xyz_movie(NumPy_Array, Z_in, file_name):
    #reorinate the array, i = atom #, x = xyz coordinate, t = time-step
    #NumPy_Array = np.einsum('ixt -> tix', NumPy_Array) 
    
    #http://www.chm.bris.ac.uk/~paulmay/temp/pcc/xyz.htm
    ff = open(file_name + ".xyz", "w")
    for index, element in enumerate(NumPy_Array):
        ff.write(str(len(element)) + "\n")
        ff.write(str(index) + "\n")
        ff.write(print2darray_tostring(element, Z_in))
        ff.write("\n")
    ff.close() 
    
    return None

def atomicname_to_Z(Z_atomicname):
    Element_Names = {'e ':0,  0:'e ',  'H ':  1,  1:'H ',  'He':  2,  2:'He',
                  'Li':  3,  3:'Li',  'Be':  4,  4:'Be',  'B ':  5,  5:'B ',  'C ':  6,  6:'C ',  'N ':  7,  7:'N ',  'O ':  8,  8:'O ',  'F ':  9,  9:'F ',  'Ne': 10, 10:'Ne', 
                  'Na': 11, 11:'Na',  'Mg': 12, 12:'Mg',  'Al': 13, 13:'Al',  'Si': 14, 14:'Si',  'P ': 15, 15:'P ',  'S ': 16, 16:'S ',  'Cl': 17, 17:'Cl',  'Ar': 18, 18:'Ar',
                  'K ': 19, 19:'K ',  'Ca': 20, 20:'Ca',  'Sc': 21, 21:'Sc',  'Ti': 22, 22:'Ti',  'V ': 23, 23:'V ',  'Cr': 24, 24:'Cr',  'Mn': 25, 25:'Mn',  'Fe': 26, 26:'Fe',  'Co': 27, 27:'Co',  'Ni': 28, 28:'Ni',  'Cu': 29, 29:'Cu',  'Zn': 30, 30:'Zn',  'Ga': 31, 31:'Ga',  'Ge': 32, 32:'Ge',  'As': 33, 33:'As',  'Se': 34, 34:'Se',  'Br': 35, 35:'Br',  'Kr': 36, 36:'Kr', 
                  'Rb': 37, 37:'Rb',  'Sr': 38, 38:'Sr',  'Y ': 39, 39:'Y ',  'Zr': 40, 40:'Zr',  'Nb': 41, 41:'Nb',  'Mo': 42, 42:'Mo',  'Tc': 43, 43:'Tc',  'Ru': 44, 44:'Ru',  'Rh': 45, 45:'Rh',  'Pd': 46, 46:'Pd',  'Ag': 47, 47:'Ag',  'Cd': 48, 48:'Cd',  'In': 49, 49:'In',  'Sn': 50, 50:'Sn',  'Sb': 51, 51:'Sb',  'Te': 52, 52:'Te',  'I ': 53, 53:'I ',  'Xe': 54, 54:'Xe',
                  'Cs': 55, 55:'Cs',  'Ba': 56, 56:'Ba',  'La': 57, 57:'La',  'Ce': 58, 58:'Ce',  'Pr': 59, 59:'Pr',  'Nd': 60, 60:'Nd',  'Pm': 61, 61:'Pm',  'Sm': 62, 62:'Sm',  'Eu': 63, 63:'Eu',  'Gd': 64, 64:'Gd',  'Tb': 65, 65:'Tb',  'Dy': 66, 66:'Dy',  'Ho': 67, 67:'Ho',  'Er': 68, 68:'Er',  'Tm': 69, 69:'Tm',  'Yb': 70, 70:'Yb',  'Lu': 71, 71:'Lu',  'Hf': 72, 72:'Hf',  'Ta': 73, 73:'Ta',  'W ': 74, 74:'W ',  'Re': 75, 75:'Re',  'Os': 76, 76:'Os',  'Ir': 77, 77:'Ir',  'Pt': 78, 78:'Pt',  'Au': 79, 79:'Au', 'Hg': 80, 80:'Hg',  'Tl': 81, 81:'Tl',  'Pb': 82, 82:'Pb',  'Bi': 83, 83:'Bi',  'Po': 84, 84:'Po', 'At': 85, 85:'At', 'Rn': 86, 86:'Rn', 
                  'Fr': 87, 87:'Fr',  'Ra': 88, 88:'Ra',  'Ac': 89, 89:'Ac',  'Th': 90, 90:'Th',  'Pa': 91, 91:'Pa',  'U ': 92, 92:'U ',  'Np': 93, 93:'Np',  'Pu': 94, 94:'Pu',  'Am': 95, 95:'Am',  'Cm': 96, 96:'Cm',  'Bk': 97, 97:'Bk',  'Cf': 98, 98:'Cf',  'Es': 99, 99:'Es',  'Fm':100,100:'Fm',  'Md':101,101:'Md',  'No':102,102:'No',  'Lr':103,103:'Lr',  'Rf':104,104:'Rf',  'Db':105,105:'Db',  'Sg':106,106:'Sg',  'Bh':107,107:'Bh',  'Hs':108,108:'Hs',  'Mt':109,109:'Mt',  'Ds':110,110:'Ds',  'Rg':111,111:'Rg', 'Cn':112,112:'Cn',  'Nh':113,113:'Nh',  'Fl':114,114:'Fl',  'Mc':115,115:'Mc',  'Lv':116,116:'Lv', 'Ts':117,117:'Ts', 'Og':118,118:'Og'}

    Z = np.zeros(len(Z_atomicname), dtype=int)
    for i in range(len(Z_atomicname)):
        Z[i] = Element_Names[(Z_atomicname[i]).strip()]

    return Z

def get_massses(Z):
    Z_mass = np.array([ 1.,   1837.,   7297.,  
                   12650., 16427.,  19705.,  21894.,  25533.,  29164.,  34631., 36785.,  
                   41908., 44305.,  49185.,  51195.,  56462.,  58441.,  64621.,  72820.,  
                   71271., 73057.,  81949.,  87256.,  92861.,  94782., 100145., 101799., 107428., 106990., 115837., 119180., 127097., 132396., 136574., 143955., 145656., 152754.,
                   155798., 159721., 162065., 166291., 169357., 174906., 176820., 184239., 187586., 193991., 196631., 204918., 209300., 216395., 221954., 232600., 231331., 239332., 
                   242270., 250331., 253208., 255415., 256859., 262937., 264318., 274089., 277013., 286649., 289702., 296219., 300649., 304894., 307947., 315441., 318945., 325367., 329848., 335119., 339434., 346768., 350390., 355616., 359048., 365656., 372561., 377702., 380947., 380983., 382806., 404681.,
                   406504., 411972., 413795., 422979., 421152., 433900., 432024., 444784., 442961., 450253., 450253., 457545., 459367., 468482., 470305., 472128., 477596., 486711., 492179., 490357., 492179., 492179., 506763., 512231., 512231., 519523., 521346., 526814., 526814., 534106., 534106., 535929.])

    return Z_mass[Z]


