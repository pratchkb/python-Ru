import numpy as np
import math


# Angle between two vectors
def angle(ax,ay,az,bx,by,bz):
    """Calculates the angle between any 2 vectors"""
    normA = np.sqrt(ax**2+ay**2+az**2)
    normB = np.sqrt(bx**2+by**2+bz**2)
    costheta = ( ax*bx + ay * by + az * bz)/(normA * normB)
    theta = math.degrees(math.acos(costheta))
    return theta


# bond vector calculate
def bond_vector(atom1_coord, atom2_coord):
    """ Calculates any bond vector"""
    x_dist = atom1_coord[0] - atom2_coord[0]
    y_dist = atom1_coord[1] - atom2_coord[1]
    z_dist = atom1_coord[2] - atom2_coord[2]  
    return (x_dist, y_dist, z_dist)


# Center of Coordinates Calculation
def COcoord(atom1_coord, atom2_coord):
    """ Calculates the center of coordinates of two atoms"""
    x_dist_c = (atom1_coord[0] + atom2_coord[0])/2
    y_dist_c = (atom1_coord[1] + atom2_coord[1])/2
    z_dist_c = (atom1_coord[2] + atom2_coord[2])/2
    return (x_dist_c, y_dist_c, z_dist_c)


# function that corrects PBCs along x-direction
def correctPBCx(atom1_coord, atom2_coord):
    """Fixes bonds stradled along x-axis.
    Need to define half_box_x, full_box_x and half_box_xcoord 
    in the main program."""
    x_dist = abs(atom1_coord[0] - atom2_coord[0])
    if x_dist > half_box_x:
        if atom2_coord[0] < half_box_xcoord:
            atom2_coord[0] = atom2_coord[0] + full_box_x
            return atom2_coord
        else:
            atom2_coord[0] = atom2_coord[0] - full_box_x
            return atom2_coord
    else:
        return atom2_coord


# function that corrects PBCs along y-direction
def correctPBCy(atom1_coord, atom2_coord):
    """Fixes bonds stradled along y-axis.
    Need to define half_box_y, full_box_y and half_box_ycoord 
    in the main program."""
    y_dist = abs(atom1_coord[1] - atom2_coord[1])
    if y_dist > half_box_y:
        if atom2_coord[1] < half_box_ycoord:
            atom2_coord[1] = atom2_coord[1] + full_box_y
            return atom2_coord
        else:
            atom2_coord[1] = atom2_coord[1] - full_box_y
            return atom2_coord
    else:
        return atom2_coord

# function that corrects PBCs along z-direction
def correctPBCz(atom1_coord, atom2_coord):
    """Fixes bonds stradled along z-axis.
    Need to define half_box_z, full_box_z and half_box_zcoord 
    in the main program."""
    z_dist = abs(atom1_coord[2] - atom2_coord[2])
    if z_dist > half_box_z:
        if atom2_coord[2] < half_box_zcoord:
            atom2_coord[2] = atom2_coord[2] + full_box_z
            return atom2_coord
        else:
            atom2_coord[2] = atom2_coord[2] - full_box_z
            return atom2_coord
    else:
        return atom2_coord


