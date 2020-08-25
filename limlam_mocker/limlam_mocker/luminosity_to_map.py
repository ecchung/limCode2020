from __future__ import absolute_import, print_function
import numpy as np
from  .tools import *
from . import debug

#@timeme
def Lco_to_Lmap(halos,map):
    """
    Bins luminosities into 3D data cube

    Parameters
    ----------
    halos : class
        Contains all halo information (position, redshift, etc..) 
    map : class
       contains all information about the map that the halos will be binned into

    Returns
    -------
    maps :
        The 3D data cube of brightness temperature
    """
    if debug.verbose: print('\n\tBinning halos into map')

    ### Calculate line freq from redshift
    halos.nu  = map.nu_rest/(halos.redshift+1)       
    
    # Transform from Luminosity to Temperature
    #halos.Tco = T_line(halos, map)

    # flip frequency bins because np.histogram needs increasing bins
    bins3D = [map.pix_binedges_x, map.pix_binedges_y, map.nu_binedges[::-1]]

    # bin in RA, DEC, NU_obs
    maps, edges = np.histogramdd( np.c_[halos.ra, halos.dec, halos.nu], 
                                  bins    = bins3D,
                                  weights = halos.Lco )
    # flip back frequency bins
    return maps[:,:,::-1]
    
def Ngal_to_map(halos,map):
    """
    Bins halo number counts into a 3d data cube

    Parameters
    ----------
    halos : class
        Contains all halo information (position, redshift, etc..) 
    map : class
       contains all information about the map that the halos will be binned into

    Returns
    -------
    maps :
        The 3D data cube of halo number counts
    """
    if debug.verbose: print('\n\tBinning halos into map')

    ### Calculate line freq from redshift
    halos.nu  = map.nu_rest/(halos.redshift+1)

    # flip frequency bins because np.histogram needs increasing bins
    bins3D = [map.pix_binedges_x, map.pix_binedges_y, map.nu_binedges[::-1]]

    # bin in RA, DEC, NU_obs
    maps, edges = np.histogramdd( np.c_[halos.ra, halos.dec, halos.nu], 
                                  bins    = bins3D)
    # flip back frequency bins
    return maps[:,:,::-1]


def T_line(halos, map):
    """
    The line Temperature in Rayleigh-Jeans limit
    T_line = c^2/2/kb/nuobs^2 * I_line

     where the Intensity I_line = L_line/4/pi/D_L^2/dnu
        D_L = D_p*(1+z), I_line units of L_sun/Mpc^2/Hz

     T_line units of [L_sun/Mpc^2/GHz] * [(km/s)^2 / (J/K) / (GHz) ^2] * 1/sr
        = [ 3.48e26 W/Mpc^2/GHz ] * [ 6.50966e21 s^2/K/kg ] 
        = 2.63083e-6 K = 2.63083 muK 
    """ 
    convfac = 2.63083
    Tco     = 1./2*convfac/halos.nu**2 * halos.Lco/4/np.pi/halos.chi**2/(1+halos.redshift)**2/map.dnu/map.Ompix

    return Tco

#@timeme
def save_maps(map):
    """
    save 3D data cube in .npz format, including map header information
    """
    if debug.verbose: print('\n\tSaving Map Data Cube to\n\t\t',map.output_file)
    np.savez(map.output_file,
             fov_x=map.fov_x, fov_y=map.fov_y,
             pix_size_x=map.pix_size_x, pix_size_y=map.pix_size_y,
             npix_x=map.npix_x, npix_y=map.npix_y,
             map_pixel_ra    = map.pix_bincents_x,
             map_pixel_dec   = map.pix_bincents_y,
             map_frequencies = map.nu_bincents,
             map_cube        = map.maps)

    return
