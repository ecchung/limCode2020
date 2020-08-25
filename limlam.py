'''
Module for using George Stein and Dongwoo Chung's limlam_mocker simulation
code in the lim structures
'''

import numpy as np
import scipy as sp
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
import h5py

from line_obs import LineObs
from _utils import cached_property,get_default_params,check_params,check_model
import _vid_tools as vt

import limlam_mocker.limlam_mocker as llm
from limlam_mocker.limlam_mocker.tools import *
from scipy.ndimage import gaussian_filter


class LimLam(LineObs):
    '''
    An object which generates limlam_mocker data cubes and computes various
    quantities from them.
    
    A subclass of LineObs, so the cosmology, emission model, and size of the
    simulated volumes are determined using the same syntax as LineModel and
    LineObs.  Outputs have astropy units added to them for ease of comparison
    with analytical results.
    
    INPUT PARAMETERS:
    catalogue_file:   String containing location of peak-patch catalogue, in
                      the form of a .npz file.  Some catalogues can be found
                      at http://cita.utoronto.ca/~gstein/data/CO/COMAP_fullvolume/peak-patch-runs/

    map_output_file:  Location to store file with output of limlam_mocker
    
    scatter_seed:     Seed for generating lognormal random scatter. If 'None',
                      scatter is fully random.  Assigning a number will
                      generate the same scatter every time.
                      
    randomize_halos   Randomize the order of the list of halo luminosities if
                      True.  Effectively assigns the same average bias to each
                      halo.
    '''
    
    #filename = 'limlam_mocker/catalogues/default_catalogue.npz'        # George's default catalogue
    filename = 'limlam_mocker/catalogues/Martines_galaxy_catalogue.h5'  # Martine's catalogue
    
    def __init__(self,
              catalogue_file=filename,
              map_output_file='limlam_output.npz',
              scatter_seed=None,
              randomize_halos=False,
              saveHalos=False,
              haloType='all',
              halos_nocutoff=False,
              **obs_kwargs):
        
        # Initialize LineObs and LineMoedl
        LineObs.__init__(self,**obs_kwargs)
        
        # Simulation parameters
        self._sim_params = locals()
        self._sim_params.pop('self')
        self._sim_params.pop('obs_kwargs')
        self._default_sim_params = get_default_params(LimLam.__init__)
        
        for key in self._sim_params:
            setattr(self,key,self._sim_params[key])
         
        # Sims currently only available for ML models
        '''
        if self.model_type!='ML':
            raise ValueError('Limlam sims only available for ML models')
        '''
        # Combine sim_params with obs_params
        self._input_params.update(self._sim_params)
        self._default_params.update(self._default_sim_params)
        
    #####################
    # Limlam parameters #
    #####################
    
    @cached_property
    def limlam_params(self):
        '''
        Translates LineModel parameter names into limlam_mocker names found
        in params.py
        '''
        
        par = empty_table()
        
        par.nmaps = self.Nch
        # Assumes square survey field
        par.fov_x = (np.sqrt(self.Omega_field).to(u.deg)).value
        par.fov_y = (np.sqrt(self.Omega_field).to(u.deg)).value
        par.npix_x = self.Nside
        par.npix_y = self.Nside
        par.nu_rest = (self.nu.to(u.GHz)).value
        par.nu_i = ((self.nuObs+self.Delta_nu/2.).to(u.GHz)).value
        par.nu_f = ((self.nuObs-self.Delta_nu/2.).to(u.GHz)).value
        par.map_output_file = self.map_output_file
        
        return par
    
    ###########################
    # Generate simulated maps #
    ###########################
    
    @cached_property
    def mapinst(self):
        '''
        Generates map object from input parameters
        '''
        return llm.params_to_mapinst(self.limlam_params)
        
    @cached_property # CLARA MADE CHANGE!!
    def halo_info(self):
        '''
        Reads in information about halos from catalogue_file
        '''     
        if self.catalogue_file.endswith('.h5'):                         #### for .h5 catalogues
            return  h5py.File(self.catalogue_file, 'r')                 ####
        else:
            return np.load(self.catalogue_file,allow_pickle=True)
        
    def check_cosmo(self,sim_cosmo,cosmo_model):
        '''
        Function to compare a limlam cosmo object to an astropy FlatLambdaCDM
        object.  Returns True if parameters match, False if any do not match
        '''
        sim_H0 = sim_cosmo.h*100*u.km/(u.s*u.Mpc)
        
        x = (sim_cosmo.Omega_M==cosmo_model.Om0,
             sim_cosmo.Omega_B==cosmo_model.Ob0,
             sim_H0==cosmo_model.H0)
        
        return all(x)
             
        
    @cached_property
    def limlam_cosmo(self):
        '''
        Cosmological model used in simulations.  If force_sim_cosmo is set to
        True and the LineObs object has a different cosmo_model, this function
        also updates the cosmo_model to match the simulations, as it is easier
        to alter the analytic model than to run new halo catalogues.
        
        Note that George's load_peakpatch_catalogue functions have been
        modified slightly from the versions on git.
        '''
        cosmo = llm.load_peakpatch_catalogue_cosmo(self.halo_info)
        
        if not self.check_cosmo(cosmo,self.h.cosmo):
            print("Input cosmological model does not match simulations")
                
        return cosmo
        
    @cached_property
    def _cosmo_flag(self):
        '''
        True if simulated cosmology does not match input to analytic model
        '''
        if not self.check_cosmo(self.limlam_cosmo,self.h.cosmo):
            return True
        else:
            return False
        
    @cached_property
    def halos(self):
        '''
        Returns all halo information (position, redshift, etc.). List of halos
        is culled to remove undesired halos.
        
        Note that George's load_peakpatch_catalogue functions have been
        modified slightly from the versions on git.
        '''
        if self.catalogue_file.endswith('.h5'):
            filetype = '.h5'
        else:
            filetype = '.npz'
        unculled_halos = llm.load_peakpatch_catalogue(self.halo_info, filetype=filetype, saveHalos=self.saveHalos)
        min_mass = self.Mmin.to(u.Msun).value
        max_mass = self.Mmax.to(u.Msun).value
        if self.halos_nocutoff:
            return unculled_halos
        else:
            return llm.cull_peakpatch_catalogue(unculled_halos,min_mass,max_mass,self.mapinst,haloType=self.haloType)
        
    @cached_property
    def L_halos(self):
        '''
        List of the line luminosities of each halo.
        
        The original github version of Mhalo_to_Lco has been modified to use
        the ML models from mass_luminosity.py
        '''
        if self.model_type=='ML':
            L =  llm.Mhalo_to_Lline(self.halos,self.model_type,
                   self.model_name,self.model_par,
                   sigma_scatter=self.sigma_scatter,
                   scatter_seed=self.scatter_seed,
                   randomize_halos=self.randomize_halos)
        elif self.model_type=='LF':
            L =  llm.Mhalo_to_Lline(self.halos,self.model_type,
                   self.model_name,self.model_par,
                   sigma_scatter=self.sigma_scatter,
                   scatter_seed=self.scatter_seed,L=self.L,dndL=self.dndL,
                   randomize_halos=self.randomize_halos)
                
        self.halos.Lco = L.to(u.Lsun).value
        return L
        
    @cached_property
    def M_halos(self):
        '''
        List of the masses of each halo
        '''
        return self.halos.M*u.Msun
    

    @cached_property
    def maps(self):
        '''
        Generate maps from computed halo luminosities
        '''
        # Make sure L_halos has been called, to ensure that halos has an
        # Lco attribute
        L_halos = self.L_halos
        
        maps = llm.Lco_to_Lmap(self.halos,self.mapinst)
        self.mapinst.maps = maps
        return maps*u.Lsun*self.XLT
    
    @cached_property
    def maps_smoothed(self):
        '''
        Map smoothed by Gaussian beam
        '''
        for ii in range(0,self.Nch):
            map = scipy.ndimage
        
    @cached_property
    def Ngal_maps(self):
        '''
        Generate map of galaxy number counts
        '''
        return llm.Ngal_to_map(self.halos,self.mapinst)
        
    @cached_property
    def noise_added_map(self): # added by Clara
        '''
        Replace 'm' with 'self' in limlam.py
        '''
        map = self.maps
        sm_map = gaussian_filter(map, [1,1,0])
        noise_sigma  = self.sigma_N/np.sqrt(self.tobs*self.Nfeeds)
        noise_map    = np.random.normal(0, noise_sigma.to(u.Jy/u.sr).value, map.shape)
        print(noise_sigma.to(u.Jy/u.sr).value)
        sm_noise_map = sm_map + noise_map
        return sm_noise_map
    
    ##################
    # Map statistics #
    ##################
    
    @cached_property
    def _Pk_dict(self):
        '''
        Dictionary containing k, Pk, and Nmodes computed from simulation. This
        function only exists to output the results of limlam_mocker in a way
        which can be easily accessed through cached_properties
        '''
        # Make sure the maps property has been called, to ensure that mapinst
        # has a map attribute
        maps = self.maps
        
        return llm.map_to_pspec(self.mapinst,self.limlam_cosmo)
        
    @cached_property
    def k_sim(self):
        '''
        Wavenumbers for simulated power spectrum.  
        '''
        return self._Pk_dict['k']/u.Mpc
    
    @cached_property
    def Pk_sim(self):
        '''
        Simulated line power spectrum
        '''
        return (self._Pk_dict['Pk']*self.XLT**2).value*self.Pk.unit
        
    @cached_property
    def cv_sim(self):
        '''
        Sample variance error on Pk from simulations
        '''
        return self.Pk_sim/np.sqrt(self._Pk_dict['nmodes'])
        
    @cached_property
    def dndM_sim(self):
        '''
        Simulated halo mass function
        '''
        Medge = vt.binctr_to_binedge_log(self.M)
        dM = np.diff(Medge)
        dNdM = np.histogram(self.M_halos,bins=Medge)[0]/dM
        return dNdM/self.Vfield
        
    @cached_property
    def nbar_sim(self):
        '''
        Mean number density of simulated halos
        '''
        return self.M_halos.size/self.Vfield
        
    @cached_property
    def dndL_sim(self):
        '''
        Simulated halo luminosity function
        '''
        Ledge = vt.binctr_to_binedge_log(self.L)
        dL = np.diff(Ledge)
        dNdL = np.histogram(self.L_halos,bins=Ledge)[0]/dL
        return dNdL/self.Vfield
        
    @cached_property
    def Tmean_sim(self):
        '''
        Mean intensity of simulated volume
        '''
        return self.maps.mean()
        
    @cached_property
    def Tmean_halos(self):
        '''
        Mean intensity estimated by summing over all halos
        '''
        return self.CLT*self.L_halos.sum()/self.Vfield
        
    @cached_property
    def Pshot_halos(self):
        '''
        Shot noise power estimated by summing over all halos
        '''
        return self.CLT**2*(self.L_halos**2).sum()/self.Vfield
        
    @cached_property
    def Bi_sim(self):
        '''
        Simulated VID
        '''
        return np.histogram(self.maps,self.Tedge_i)[0]
        
    @cached_property
    def PofN_sim(self):
        '''
        Computes number count PDF from simulated volume.  To be compared
        to self.PofN
        '''
        # Number count bin edges
        Ngal_edge = np.linspace(-0.5,self.Ngal_max+0.5,self.Ngal_max+2)
        return np.histogram(self.Ngal_maps,Ngal_edge)[0]/self.Nvox
        
        
        
        
###############################
# Set cosmology to match sims #
###############################
def set_sim_cosmo(limlam):
    '''
    Takes a LIMLAM object and updates it so that the cosmologies used in the
    analytic and numerical calculations match
    '''
    
    cosmo = limlam.limlam_cosmo
    new_cosmo = FlatLambdaCDM(H0=cosmo.h*100*u.km/(u.Mpc*u.s),
                              Om0=cosmo.Omega_M,Tcmb0=2.725*u.K,
                              m_nu=[0,0,0.06]*u.eV,Ob0=cosmo.Omega_B)
    limlam.update(cosmo_model=new_cosmo)
    return new_cosmo

        
