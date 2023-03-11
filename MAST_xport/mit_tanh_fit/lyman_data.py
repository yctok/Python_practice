'''
Methods to load and plot C-Mod Ly-alpha data.

sciortino, August 2020
'''
import matplotlib.pyplot as plt
plt.ion()
import numpy as np
import xarray
from scipy.interpolate import interp1d, interp2d
#sys.path.append('/home/millerma/OMFIT-source/omfit')
from omfit_classes import omfit_eqdsk, omfit_mds
import shutil, os, scipy, copy
from IPython import embed
import MDSplus
from omfit_classes.omfit_mds import OMFITmdsValue

from scipy.constants import Boltzmann as kB, e as q_electron
from scipy.optimize import curve_fit
import aurora
import sys
sys.path.append('/home/sciortino/usr/python3modules/profiletools3')
sys.path.append('/home/sciortino/usr/python3modules/eqtools3')
#sys.path.append('/home/millerma/python3modules/eqtools3')
#sys.path.append('/home/millerma/python3modules/profiletools3')
import profiletools
import eqtools

# from mitlya repo
import mtanh_fitting
import fit_2D
import tomographic_inversion as tomo


def get_cmod_kin_profs(shot,tmin,tmax, geqdsk=None, pre_shift_TS=False, force_to_zero=False, flat_err=True,
                       probes=['A'], osborne_fit=False, apply_final_sep_stretch=False):
    '''Function to load and fit modified-tanh functions to C-Mod ne and Te.

    This function is designed to be robust for operation within the construction of
    the C-Mod Ly-a database. It makes use of the profiletools package just to conveniently
    collect data and time-average it. An updated profiletools (python 3+) is available 
    from one of the github repos of sciortinof.

    Parameters
    ----------
    shot : int
        CMOD shot number
    tmin, tmax : floats
        Times in seconds for window interval of interest
    geqdsk : dict
        Dictionary containing processed EFIT geqdsk file
    pre_shift_TS : bool, opt
        If True, Thomson Scattering (TS) data are shifted based on the 2-point model after being fitted 
        by themselves. This is recommended only if TS has good coverage inside and outside of the LCFS
        and the LCFS is therefore reasonably covered by experimental data.
    force_to_zero : bool, opt
        If True, add fictitious points far into the SOL to enforce that both ne and Te must go towards
        zero at the wall.
    probes : list of str
        List of strings indicating data from which probes should be loaded if possible. 
        Possible options are 'A' (ASP probe) and 'F' (FSP probe).
    osborne_fit : bool
        If True, use the Osborne mtanh fit, otherwise use the version modified by FS.
    apply_final_sep_stretch : bool
        If True, final kinetic profiles are stretched such that Te_sep from the 2-point model is matched.

    Returns
    -------
    roa_kp : 1D array
        r/a grid
    ne : 1D array
        Electron density in units of :math:`10^{20} m^{-3}`
    ne_std : 1D array
        Uncertainties on electron density in units of :math:`10^{20} m^{-3}`
    Te : 1D array
        Electron temperature in units of :math:`keV`
    Te_std : 1D array
        Uncertainties on electron temperature in units of :math:`keV`
    p_ne : profiletools object
        Electron density object containing experimental data from all loaded diagnostics.
    p_Te : profiletools object
        Electron temperature object containing experimental data from all loaded diagnostics.
    '''
    
    if geqdsk is None:
        # the geqdsk is used only when using the Aurora radial coordinate transformations
        geqdsk = get_geqdsk_cmod(
            shot, (tmin+tmax)/2.*1e3, gfiles_loc = '/home/sciortino/EFIT/lya_gfiles/')

    # number of iterations for osborne fit to converge
    maxfev = 2000

    # identify mode to decide how to fit profile - empirical observation
    mode = identify_mode(int(shot))
    ped_width = 0.04 if mode == 'H' else 0.02 # will later be used to fit

    try:
        # require edge Thomson to be available
        p_Te= profiletools.Te(int(shot), include=['ETS'], abscissa='psinorm',t_min=tmin,t_max=tmax, remove_zeros=True)
        p_ne= profiletools.ne(int(shot), include=['ETS'], abscissa='psinorm',t_min=tmin,t_max=tmax, remove_zeros=True)
    except MDSplus.TreeNODATA:
        raise ValueError('No edge Thomson data!')
  
    try:
        equal_R = p_ne.X[:,1] == p_Te.X[:,1]
        assert np.sum(equal_R) == len(p_ne.X[:,1])
    except:
        raise ValueError('Edge Thomson rhobase differs between ne and Te')
    
    
    # if shot==1080416025:
    #     # add data from sister shot 1101014029
    #     p_Te_2= profiletools.Te(1101014029, include=['ETS'], abscissa='r/a',t_min=tmin,t_max=tmax)
    #     p_ne_2= profiletools.ne(1101014029, include=['ETS'], abscissa='r/a',t_min=tmin,t_max=tmax)
    #     p_Te.add_profile(p_Te_2)
    #     #p_Te = p_Te_2 # test sub
    #     p_ne.add_profile(p_ne_2)
    
    try:
        # try to add core Thomson, not strictly necessary 
        p_Te_CTS = profiletools.Te(int(shot), include=['CTS'],
                                   abscissa='psinorm',t_min=tmin,t_max=tmax)
        p_ne_CTS = profiletools.ne(int(shot), include=['CTS'],
                                   abscissa='psinorm',t_min=tmin,t_max=tmax)
        #p_Te.remove_points(p_Te.X[:,1]>0.9
        p_Te.add_profile(p_Te_CTS)
        p_ne.add_profile(p_ne_CTS)
    except Exception:
        pass
   
    #try:
    #    # try to add GPC and GPC2 for Te, not strictly necessary
    #    p_Te_GPC= profiletools.Te(int(shot), include=['GPC','GPC2'],
    #                              abscissa='r/a',t_min=tmin,t_max=tmax)

    #    # downsample to fewer points in time interval - make it comparable to Thomson
    #    p_Te_GPC.keep_slices(0, np.linspace(tmin,tmax, p_Te.X.shape[0]))
    #    p_Te_GPC.remove_points(np.abs(p_Te_GPC.X[:,1])>0.85)  # doubtful about ECE opacity
    #    p_Te.add_profile(p_Te_GPC)
    #except:
    #    pass


    # consider only flux surface on which points were measured, regardless of LFS or HFS
    p_Te.X=np.abs(p_Te.X)
    p_ne.X=np.abs(p_ne.X)

    # set some minimum uncertainties. Recall that units in objects are 1e20m^{-3} and keV
    p_ne.y[p_ne.y<=0.] = 0.01  # 10^18 m^-3
    p_Te.y[p_Te.y<=0.01] = 0.01 # 10 eV
    p_ne.err_y[p_ne.err_y<=0.01] = 0.01 # 10^18 m^-3
    p_Te.err_y[p_Te.err_y<=0.02] = 0.02 # 20 eV
    
    # points in the pedestal that have x uncertainties larger than 0.1 don't help at all
    # do this filtering here because filtering of err_X only works before time-averaging
    #p_ne.remove_points(np.logical_and(p_ne.X[:,1]>0.9, p_ne.err_X[:,1]>0.1))
    #p_Te.remove_points(np.logical_and(p_Te.X[:,1]>0.9, p_Te.err_X[:,1]>0.1))
    
    # time average now, before trying to add time-independent probe data
    #p_ne.time_average(weighted=True)
    #p_Te.time_average(weighted=True)
   
    p_ne.drop_axis(0)
    p_Te.drop_axis(0)

    # cleanup of low Te values
    p_Te.remove_points(np.logical_and(p_Te.X[:,0]<1.03, p_Te.y<0.015))  # TS Te should be >15 eV inside near SOL

    if shot==1080416025:
        #p_Te.err_y[p_Te.err_y<0.03] = 0.03
        #p_Te.err_y[-2:] /= 50.
        #p_ne.err_y[-2:] /= 50.

        pass

        # temporary!!!!
        #p_Te.remove_points(p_Te.X[:,0]>1)
        #p_ne.remove_points(p_ne.X[:,0]>1)
        
    #p_ne.err_y[-2:] /= 3.
    
    # use two point model to get T_sep
    Te_sep_eV = fit_2D.Teu_2pt_model(shot,tmin,tmax, geqdsk, pressure_opt = 3)
    print('Te LCFS eV', Te_sep_eV)

    
    pre_shift_TS = True
    roa_kp = np.linspace(0.0, 1.16, 300) # can't extend too much past 1.16
    
    ne_X_before = p_ne.X[:,0]
    Te_X_before = p_Te.X[:,0]

    reg = [4,8]
    if pre_shift_TS:

        p_ne, p_Te = prefit_filter(p_ne, p_Te)
    
        if force_to_zero:
            # force fits to go down in the far SOL (r/a=1.2) to make the routine more robust
            p_Te.add_data( np.array([[1.2]]), np.array([10e-3]), err_X=np.array([0.001]), err_y=np.array([0.001]))
            p_ne.add_data( np.array([[1.2]]), np.array([0.1]), err_X=np.array([0.001]), err_y=np.array([0.001]))

        # need to fit TS first 
        #min_TS_X = p_Te.X.min() if p_Te.X.min() < p_ne.X.min() else p_ne.X.min()
        #max_TS_X = p_Te.X.max() if p_Te.X.max() > p_ne.X.max() else p_ne.X.max()
        #X_TS_fit = np.linspace(min_TS_X,max_TS_X,100)
        
        idxs_ne = np.argsort(p_ne.X[:,0])
        idxs_Te = np.argsort(p_Te.X[:,0])

        mask_ne = np.full(p_ne.X[:,0].shape, True)
        mask_Te = np.full(p_Te.X[:,0].shape, True)
        
        # fit the profiles!        
        maxfev = 2000 # number of iterations for osbourne fit to converge
    
        ne, ne_popt, ne_perr, ne_chisqr = mtanh_fitting.best_osbourne(p_ne.X[idxs_ne,0][mask_ne], p_ne.y[idxs_ne][mask_ne], vals_unc=p_ne.err_y[idxs_ne][mask_ne], x_out=roa_kp, maxfev=maxfev, reg=reg, ped_width=ped_width, ne=True)
        Te, Te_popt, Te_perr, Te_chisqr = mtanh_fitting.best_osbourne(p_Te.X[idxs_Te,0][mask_Te], p_Te.y[idxs_Te][mask_Te], vals_unc=p_Te.err_y[idxs_Te][mask_Te], x_out=roa_kp, maxfev=maxfev, reg=reg, ped_width=ped_width, ne=False)
        #Te, Te_popt, Te_perr, Te_chisqr = mtanh_fitting.super_fit(p_Te.X[idxs_Te,0][mask_Te], p_Te.y[idxs_Te][mask_Te], p_Te.err_y[idxs_Te][mask_Te], x_out=roa_kp, edge_focus=True, plot=False, bounds=None)

        # apply filter depending on initial fit
        p_ne, p_Te = postfit_filter(p_ne, p_Te, ne, Te, roa_kp)

        if force_to_zero:
            # add xtra points in the far SOL again to force fits to go down
            p_Te.add_data( np.array([[1.2]]), np.array([10e-3]), err_X=np.array([0.001]), err_y=np.array([0.001]))
            p_ne.add_data( np.array([[1.2]]), np.array([0.1]), err_X=np.array([0.001]), err_y=np.array([0.001]))

        # fit them again!
        idxs_ne = np.argsort(p_ne.X[:,0])
        idxs_Te = np.argsort(p_Te.X[:,0])

        mask_ne = np.full(p_ne.X[:,0].shape, True)
        mask_Te = np.full(p_Te.X[:,0].shape, True)
        
        ne, ne_popt, ne_perr, ne_chisqr = mtanh_fitting.best_osbourne(p_ne.X[idxs_ne,0][mask_ne], p_ne.y[idxs_ne][mask_ne], vals_unc=p_ne.err_y[idxs_ne][mask_ne], x_out=roa_kp, maxfev=maxfev, reg=reg, ped_width=ped_width, ne=True)
        Te, Te_popt, Te_perr, Te_chisqr = mtanh_fitting.best_osbourne(p_Te.X[idxs_Te,0][mask_Te], p_Te.y[idxs_Te][mask_Te], vals_unc=p_Te.err_y[idxs_Te][mask_Te], x_out=roa_kp, maxfev=maxfev, reg=reg, ped_width=ped_width, ne=False)

        ## try to shift TS profiles using T_lcfs
        #xSep_TS = fit_2D.shift_profs([1],X_TS_fit,Te_TS_fit[None,:]*1e3,Te_LCFS=Te_sep_eV)
        xSep_TS = fit_2D.shift_profs([1],roa_kp,Te[None,:]*1e3,Te_LCFS=Te_sep_eV)
        print('Shifting raw TS points by r/a = {:.3f}'.format(float(1 - xSep_TS)))
        p_ne.X += 1 - xSep_TS
        p_Te.X += 1 - xSep_TS
    
        # before final fit, initalize electron pressure profile and fit that as well
        # initialize pressure profile tool object and propagate uncertainties from ne, Te
        p_pe = create_pe(p_ne, p_Te)  
        mask_pe = np.full(p_pe.X[:,0].shape, True)
        idxs_pe = np.argsort(p_pe.X[:,0])

        # now fit!
        ne, ne_popt, ne_perr, ne_chisqr = mtanh_fitting.best_osbourne(p_ne.X[idxs_ne,0][mask_ne], p_ne.y[idxs_ne][mask_ne], vals_unc=p_ne.err_y[idxs_ne][mask_ne], x_out=roa_kp, maxfev=maxfev, reg=reg, ped_width=ped_width, ne=True)
        Te, Te_popt, Te_perr, Te_chisqr = mtanh_fitting.best_osbourne(p_Te.X[idxs_Te,0][mask_Te], p_Te.y[idxs_Te][mask_Te], vals_unc=p_Te.err_y[idxs_Te][mask_Te], x_out=roa_kp, maxfev=maxfev, reg=reg, ped_width=ped_width, ne=False)
        pe, pe_popt, pe_perr, pe_chisqr = mtanh_fitting.best_osbourne(p_pe.X[idxs_pe,0][mask_pe], p_pe.y[idxs_pe][mask_pe], vals_unc=p_pe.err_y[idxs_pe][mask_pe], x_out=roa_kp, maxfev=maxfev, reg=reg, ped_width=ped_width, ne=False)

    # attempt to fetch ASP and FSP data if available

    use_pressure = False # if set to True will not shift probes and will try to match TS + SP pressure at separatrix
    
    p_ne_p, p_ne_T = None, None

    try:
        p_ne_p, p_Te_p  = fetch_edge_probes(shot, (tmin+tmax)/2., Te_sep_eV, geqdsk = geqdsk,
                                        roa_min=0.9, # 0.96,  # fetch all data before radially shifting
                                        roa_max=1.1,  # don't trust data too far into the SOL...
                                        probes=probes, shift_probes = not use_pressure) # if use pressure don't shift probes
    except:
        print('No good probe data for this shot!')

    # combine profiles and refit + shift if there is probe data
    if p_ne_p is not None: 

        if use_pressure:
            _out = match_pressure(p_ne, p_Te, p_ne_p, p_Te_p, plot=True)
            p_ne_p, p_Te_p = _out

        filt_TS = False
        filt_SP = False

        # option to filter some data above certain range
        
        if filt_TS:
            _out = filter_TS(p_ne, p_Te, p_ne_p, p_Te_p, 0.02) # cutoff in keV
            p_ne, p_Te, p_ne_p, p_Te_p = _out

        if filt_SP:
            _out = filter_SP(p_ne, p_Te, p_ne_p, p_Te_p, 0.05) # cutoff in keV
            p_ne, p_Te, p_ne_p, p_Te_p = _out
    
        # add cleaned profiles
        p_ne.add_profile(p_ne_p)
        p_Te.add_profile(p_Te_p)

        num_ne_SP = len(p_ne_p.X)
        num_Te_SP = len(p_Te_p.X)

        ne_X_before = p_ne.X[:,0]
        Te_X_before = p_Te.X[:,0]

        #if shot==1080416025:
            # remove points that are obviously bad
        #    p_Te.remove_points(np.logical_and(p_Te.X[:,0]>0.8, p_Te.y>1.25))
   
        # perform a re-prefilter - probably only does a filter of SP at this point     
        p_ne, p_Te = prefit_filter(p_ne, p_Te, TS=False)
    
        #if force_to_zero:
            # force fits to go down in the far SOL (r/a=1.2) to make the routine more robust
        #    p_Te.add_data( np.array([[1.2]]), np.array([10e-3]), err_X=np.array([0.001]), err_y=np.array([0.001]))
        #    p_ne.add_data( np.array([[1.2]]), np.array([0.1]), err_X=np.array([0.001]), err_y=np.array([0.001]))

        #lim_dom = 0 

        # Now fit:    
        idxs_ne = np.argsort(p_ne.X[:,0])
        idxs_Te = np.argsort(p_Te.X[:,0])

        #mask_ne = np.logical_and(p_ne.X[:,0][idxs_ne] > lim_dom, p_ne.X[:,0][idxs_ne] < 1.05) 
        #mask_Te = np.logical_and(p_Te.X[:,0][idxs_Te] > lim_dom, p_Te.X[:,0][idxs_Te] < 1.05) 

        mask_ne = np.full(p_ne.X[:,0].shape, True)
        mask_Te = np.full(p_Te.X[:,0].shape, True)

        # note: reg is only used in osbourne fit, not in superfit

        # fit the profiles! (again)
        maxfev = 2000 # number of iterations for osbourne fit to converge
    
        ne, ne_popt, ne_perr, ne_chisqr = mtanh_fitting.best_osbourne(p_ne.X[idxs_ne,0][mask_ne], p_ne.y[idxs_ne][mask_ne], vals_unc=p_ne.err_y[idxs_ne][mask_ne], x_out=roa_kp, maxfev=maxfev, reg=reg, ped_width=ped_width, ne=True)
        Te, Te_popt, Te_perr, Te_chisqr = mtanh_fitting.best_osbourne(p_Te.X[idxs_Te,0][mask_Te], p_Te.y[idxs_Te][mask_Te], vals_unc=p_Te.err_y[idxs_Te][mask_Te], x_out=roa_kp, maxfev=maxfev, reg=reg, ped_width=ped_width, ne=False)
    
        p_ne, p_Te = postfit_filter(p_ne, p_Te, ne, Te, roa_kp)

        #if force_to_zero:
        #    # add xtra points in the far SOL again to force fits to go down
        #    p_Te.add_data( np.array([[1.2]]), np.array([10e-3]), err_X=np.array([0.001]), err_y=np.array([0.001]))
        #    p_ne.add_data( np.array([[1.2]]), np.array([0.1]), err_X=np.array([0.001]), err_y=np.array([0.001]))
    
        # Fit again:
        idxs_ne = np.argsort(p_ne.X[:,0])
        idxs_Te = np.argsort(p_Te.X[:,0])

        #lim_dom = 0

        #mask_ne = np.logical_and(p_ne.X[:,0][idxs_ne] > lim_dom, p_ne.X[:,0][idxs_ne] < 1.05) 
        #mask_Te = np.logical_and(p_Te.X[:,0][idxs_Te] > lim_dom, p_Te.X[:,0][idxs_Te] < 1.05) 
    
        mask_ne = np.full(p_ne.X[:,0].shape, True)
        mask_Te = np.full(p_Te.X[:,0].shape, True)

        # before final fit, initalize electron pressure profile and fit that as well
        # initialize pressure profile tool object and propagate uncertainties from ne, Te
        p_pe = create_pe(p_ne, p_Te)  
        mask_pe = np.full(p_pe.X[:,0].shape, True)
        idxs_pe = np.argsort(p_pe.X[:,0])
    
        # fit them again! (again)
        ne, ne_popt, ne_perr, ne_chisqr = mtanh_fitting.best_osbourne(p_ne.X[idxs_ne,0][mask_ne], p_ne.y[idxs_ne][mask_ne], vals_unc=p_ne.err_y[idxs_ne][mask_ne], x_out=roa_kp, maxfev=maxfev, reg=reg, ped_width=ped_width, ne=True)
        Te, Te_popt, Te_perr, Te_chisqr = mtanh_fitting.best_osbourne(p_Te.X[idxs_Te,0][mask_Te], p_Te.y[idxs_Te][mask_Te], vals_unc=p_Te.err_y[idxs_Te][mask_Te], x_out=roa_kp, maxfev=maxfev, reg=reg, ped_width=ped_width, ne=False)
        pe, pe_popt, pe_perr, pe_chisqr = mtanh_fitting.best_osbourne(p_pe.X[idxs_pe,0][mask_pe], p_pe.y[idxs_pe][mask_pe], vals_unc=p_pe.err_y[idxs_pe][mask_pe], x_out=roa_kp, maxfev=maxfev, reg=reg, ped_width=ped_width, ne=False)

        apply_final_sep_shift = True
        if apply_final_sep_shift:
            
            # shift combined profiles once more
    
            xSep_all = fit_2D.shift_profs([1],roa_kp,Te[None,:]*1e3,Te_LCFS=Te_sep_eV)
            print('Combined profile shift: r/a = {:.3f}'.format(float(1 - xSep_all)))

            p_ne.X += 1 - xSep_all
            p_Te.X += 1 - xSep_all
        
            # before final fit, initalize electron pressure profile and fit that as well
            # initialize pressure profile tool object and propagate uncertainties from ne, Te
            p_pe = create_pe(p_ne, p_Te)  
            mask_pe = np.full(p_pe.X[:,0].shape, True)
            idxs_pe = np.argsort(p_pe.X[:,0])

            # now refit to get final profiles
            ne, ne_popt, ne_perr, ne_chisqr = mtanh_fitting.best_osbourne(p_ne.X[idxs_ne,0][mask_ne], p_ne.y[idxs_ne][mask_ne], vals_unc=p_ne.err_y[idxs_ne][mask_ne], x_out=roa_kp, maxfev=maxfev, reg=reg, ped_width=ped_width, ne=True)
            Te, Te_popt, Te_perr, Te_chisqr = mtanh_fitting.best_osbourne(p_Te.X[idxs_Te,0][mask_Te], p_Te.y[idxs_Te][mask_Te], vals_unc=p_Te.err_y[idxs_Te][mask_Te], x_out=roa_kp, maxfev=maxfev, reg=reg, ped_width=ped_width, ne=False)
            pe, pe_popt, pe_perr, pe_chisqr = mtanh_fitting.best_osbourne(p_pe.X[idxs_pe,0][mask_pe], p_pe.y[idxs_pe][mask_pe], vals_unc=p_pe.err_y[idxs_pe][mask_pe], x_out=roa_kp, maxfev=maxfev, reg=reg, ped_width=ped_width, ne=False)
            
    
    else:
        num_ne_SP = 0
        num_Te_SP = 0

    
    # now that functions are fit, get gradients

    grad_ne = get_fit_gradient(ne, ne_popt, roa_kp, 'osborne', grad_type='analytic', reg=reg, plot=False)
    grad_Te = get_fit_gradient(Te, Te_popt, roa_kp, 'osborne', grad_type='analytic', reg=reg, plot=False)
    grad_pe = get_fit_gradient(pe, pe_popt, roa_kp, 'osborne', grad_type='analytic', reg=reg, plot=False)

    f_ne = kinetic_profile(roa_kp, ne, grad_ne)
    f_Te = kinetic_profile(roa_kp, Te, grad_Te)
    f_pe = kinetic_profile(roa_kp, pe, grad_pe)

    num_mc = 150

    if flat_err:
        ne_std = ne*0.2 # 20% -- tanh fitting not set up to provide good uncertainties
        Te_std = Te*0.2 # 20%
        pe_std = np.sqrt(ne_std**2 + Te_std**2) # 20%
        
        grad_ne_std = grad_ne*0.2 # 20% -- tanh fitting not set up to provide good uncertainties
        grad_Te_std = grad_Te*0.2 # 20%
        grad_pe_std = np.sqrt(grad_ne_std**2 + grad_Te_std**2) # 20%

        # if did not do MC iteration, just set the mean coefficients as the fit coefficients
        ne_mean_popt = ne_popt
        Te_mean_popt = Te_popt
        pe_mean_popt = pe_popt
        

    else: # this is getting a bit messy so might just put some of these into an object
        f_ne.perturb_mc_err(p_ne, num_mc, maxfev, ped_width, reg, ne=True, verbose=False, plot=False)
        f_Te.perturb_mc_err(p_Te, num_mc, maxfev, ped_width, reg, ne=False, verbose=False, plot=False)
        f_pe.perturb_mc_err(p_pe, num_mc, maxfev, ped_width, reg, ne=False, verbose=False, plot=False)
            

#        _out = perturb_mc_err(num_mc, p_ne, p_Te, p_pe,  ne, Te, pe, roa_kp, maxfev, ped_width, reg=reg, remove_outliers=False, verbose=False, plot=True)
#        ne_mean, ne_std, ne_MC_fits, Te_mean, Te_std, Te_MC_fits, pe_mean, pe_std, pe_MC_fits, \
#            grad_ne_mean, grad_ne_std, grad_Te_mean, grad_Te_std, grad_pe_mean, grad_pe_std = _out

#        ne_mean, ne_std, Te_mean, Te_std, pe_mean, pe_std, \
#            grad_ne_mean, grad_ne_std, grad_Te_mean, grad_Te_std, grad_pe_mean, grad_pe_std, = \
#            ne_mean[0], ne_std[0], Te_mean[0], Te_std[0], pe_mean[0], pe_std[0], \
#            grad_ne_mean[0], grad_ne_std[0], grad_Te_mean[0], grad_Te_std[0], grad_pe_mean[0], grad_pe_std[0]

        # now refit the mean and save coefficients of mean
    
        ne_mean_fit, ne_mean_popt, ne_mean_perr, ne_mean_chisqr = mtanh_fitting.best_osbourne(roa_kp, f_ne.fit_mean, vals_unc=f_ne.fit_std, x_out=f_ne.x, maxfev=maxfev, reg=reg, ped_width=ped_width, ne=True)
        Te_mean_fit, Te_mean_popt, Te_mean_perr, Te_mean_chisqr = mtanh_fitting.best_osbourne(roa_kp, f_Te.fit_mean, vals_unc=f_Te.fit_std, x_out=f_Te.x, maxfev=maxfev, reg=reg, ped_width=ped_width, ne=False)
        pe_mean_fit, pe_mean_popt, pe_mean_perr, pe_mean_chisqr = mtanh_fitting.best_osbourne(roa_kp, f_pe.fit_mean, vals_unc=f_pe.fit_std, x_out=f_pe.x, maxfev=maxfev, reg=reg, ped_width=ped_width, ne=False)
    
        plot_fit = False
        if plot_fit:
            fig, ax = plt.subplots(2)
            ax[0].plot(p_ne.X[:,0], p_ne.y, 'o')
            ax[0].plot(f_ne.x, f_ne.y)
            ax[0].plot(f_ne.x, f_ne.fit_mean)

            ax[1].plot(p_Te.X[:,0], p_Te.y, 'o')
            ax[1].plot(f_Te.x, f_Te.y)
            ax[1].plot(f_Te.x, f_Te.fit_mean)

            ax[0].legend(['raw','fit','mean_fit'])
    
    # store fit coefficients
  
    f_ne.save_coefs(ne_popt, ne_mean_popt)
    f_Te.save_coefs(Te_popt, Te_mean_popt)
    f_pe.save_coefs(pe_popt, pe_mean_popt)
 

    '''
    iter_range = np.arange(500)+1
    _out = perturb_mc_err(iter_range, p_ne, p_Te, ne, Te, roa_kp, idxs_ne, idxs_Te, mask_ne, mask_Te, maxfev, remove_outliers=True, verbose=True, plot=False)
    ne_err_range, Te_err_range = _out

    fig, ax = plt.subplots(2, sharex=True)
    ax[0].plot(iter_range, ne_err_range[:,49])
    ax[0].plot(iter_range, ne_err_range[:,99])
    ax[0].plot(iter_range, ne_err_range[:,149])
    ax[0].plot(iter_range, ne_err_range[:,199])
    ax[0].plot(iter_range, ne_err_range[:,249])
    ax[0].plot(iter_range, ne_err_range[:,299])
    
    ax[1].plot(iter_range, Te_err_range[:,49]*1e3)
    ax[1].plot(iter_range, Te_err_range[:,99]*1e3)
    ax[1].plot(iter_range, Te_err_range[:,149]*1e3)
    ax[1].plot(iter_range, Te_err_range[:,199]*1e3)
    ax[1].plot(iter_range, Te_err_range[:,249]*1e3)
    ax[1].plot(iter_range, Te_err_range[:,299]*1e3)
   
    ax[0].set_xlabel('Number of MC iterations')
    ax[0].set_xlabel('$\Delta n_{e} [10^{20} m^{-3}]$')
    ax[1].set_xlabel('$\Delta T_{e} [eV]$')

    ax[0].legend(['r/a = 0.19','r/a = 0.38','r/a = 0.58','r/a = 0.77', 'r/a = 0.97','r/a = 1.16'])
    
    fig.suptitle(f'C-Mod shot {shot}')
    '''

    p_ne.remove_points(p_ne.X[:,0]>1.1) # remove artificial point
    p_Te.remove_points(p_Te.X[:,0]>1.1) # remove artificial point
    
    ne_X_after = p_ne.X[:,0]
    Te_X_after = p_Te.X[:,0]
    
    ne_pts_removed = np.setdiff1d(ne_X_before,ne_X_after)
    Te_pts_removed = np.setdiff1d(Te_X_before,Te_X_after)
   
    if p_ne_p is not None:
        for pt in ne_pts_removed:
            if pt in p_ne_p.X:
                num_ne_SP -= 1
        for pt in Te_pts_removed:
            if pt in p_Te_p.X:
                num_Te_SP -= 1

    nu_star = calc_nu_star(ne, Te)

    plot_fit = False
    if plot_fit:
        p_ne.plot_data()
        #plt.xlim([0.86,1.16])
        #plt.ylim([-0.01, 0.75])
        #plt.gca().plot(roa_kp, ne2)
        plt.gca().plot(roa_kp, ne)
        #plt.legend(['super_fit', 'super_fit_osborne'])

        p_Te.plot_data()
        #plt.xlim([0.86,1.16])
        #plt.ylim([-0.01, 0.5])
        #plt.gca().plot(roa_kp, Te2)
        plt.gca().plot(roa_kp, Te)
        #plt.legend(['super_fit', 'super_fit_osborne'])

    num_SP = {}
    num_SP['ne'] = num_ne_SP; num_SP['Te'] = num_Te_SP

    # plot preshift and postshift
    #fig,ax = plt.subplots(2, sharex=True)
    
    #preshift
    #plot_pre = ax[0].errorbar(p_ne.X[:,0], p_ne.y, p_ne.err_y, fmt='o')
    #color_pre = plot_pre[0].get_color()
    #ax[1].errorbar(p_Te.X[:,0], p_Te.y, p_Te.err_y, fmt='o')

    #plt.show()
       
    # check if there exist any data point in the SOL
    if np.all(p_ne.X[:,0]<0.99):
        #raise ValueError(f'No ne data points  for r/a<0.99 in shot {shot}!')
        print(f'No SOL ne data points for r/a>0.99 in shot {shot}!')
    if np.all(p_Te.X[:,0]<0.99):
     
        #raise ValueError(f'No SOL Te data points for r/a<0.99 in shot {shot}!')
        print(f'No SOL Te data points for r/a>0.99 in shot {shot}!')

    #posthift
    #plot_post = ax[0].errorbar(p_ne.X[:,0], p_ne.y, p_ne.err_y, color='blue', fmt='o')
    #color_post = plot_post[0].get_color()
    #ax[0].legend(['pre-shift', 'post-shift'])

    #ax[1].errorbar(p_Te.X[:,0], p_Te.y, p_Te.err_y, color='red', fmt='o')
    #ax[1].legend(['pre-shift', 'post-shift'])

    #ax[0].plot(roa_kp - (1-roa_of_TeSep), ne, color=color_pre)
    #ax[0].plot(roa_kp, ne, color='blue')
    #ax[1].plot(roa_kp - (1-roa_of_TeSep), Te, color=color_pre)
    #ax[1].plot(roa_kp, Te, color='red')
   
    #ax[0].axvline(roa_of_TeSep, linestyle='--', color=color_pre)
    #ax[1].axvline(roa_of_TeSep, linestyle='--', color=color_pre)
    #ax[0].axvline(1, linestyle='--', color=color_post)
    #ax[1].axvline(1, linestyle='--', color=color_post)

    #ax[0].set_ylabel(r'$n_e$ [$cm^{-3}$]', fontsize=14)
    #ax[0].tick_params(axis='y', labelsize=14)

    #ax[1].set_ylabel(r'$T_e$ [$keV$]', fontsize=14)
    #ax[1].set_xlabel(r'$r/a$', fontsize=14)
    #ax[1].tick_params(axis='x', labelsize=14)
    #ax[1].tick_params(axis='y', labelsize=14)

    #ax[0].set_xlim([0.86,1.13])
    #ax[0].set_ylim([-0.05, 1.75])
    #ax[1].set_ylim([-0.05, 1.1])

    #plt.show()
    
    # before returning things, recreate the pe profile
    p_pe = create_pe(p_ne, p_Te)  
    
    # output fits + profiletool objects for ne and Te so that experimental data points are passed too
    return f_ne, f_Te, f_pe, p_ne, p_Te, p_pe, num_SP 

            

def prefit_filter(p_ne, p_Te, TS=True):

    # apply these minimum conditions to all TS data (probes sometimes have smaller error bars)
    if TS:
        p_ne.y[p_ne.y<=0.] = 0.01  # 10^18 m^-3
        p_Te.y[p_Te.y<=0.01] = 0.01 # 10 eV
        p_ne.err_y[p_ne.err_y<=0.1] = 0.1 # 10^19 m^-3
        p_Te.err_y[p_Te.err_y<=0.02] = 0.02 # 20 eV

    # remove points with excessively large error bars
    p_ne.remove_points(p_ne.err_y>1) # 10^20 m^-3
    p_Te.remove_points(np.logical_and(p_Te.X[:,0]>0.98, p_Te.err_y>0.2)) # max 200 eV of uncertainty in the pedestal
    
    # Remove points with too high values in the SOL:
    p_ne.remove_points(np.logical_and(p_ne.X[:,0]>1.0, p_ne.y>1.5)) # 5e20
    p_Te.remove_points(np.logical_and(p_Te.X[:,0]>1.0, p_Te.y>0.25)) # 250 eV

    # Remove ne points with too high uncertainty in the pedestal and SOL:
    p_ne.remove_points(np.logical_and(p_ne.X[:,0]>0.9, p_ne.err_y>0.3)) # 3e19 m^-3
    #if shot==1100308004:
    #    p_ne.remove_points(np.logical_and(p_ne.X[:,0]>1.0, p_ne.err_y>0.2)) # 2e19 m^-3 , less in the SOL

    # Remove Te points with too high uncertainty in the SOL (better not to filter in the pedestal, high variability)
    p_Te.remove_points(np.logical_and(p_Te.X[:,0]>1.0, p_Te.err_y>0.1))  # 100 eV

    # trivial clean up of Te in the pedestal
    p_Te.remove_points(np.logical_and(p_Te.X[:,0]>0.9, p_Te.err_y>0.5))  # 500 eV
    p_Te.remove_points(np.logical_and(p_Te.X[:,0]<0.98, p_Te.y<0.05))  # Te inside r/a=0.98 must always be >50 eV
    p_Te.remove_points(np.logical_and(p_Te.X[:,0]<0.95, p_Te.y<0.1))  # Te inside r/a=0.95 must always be >100 eV
   
    # clean up bad core points
    p_Te.remove_points(np.logical_and(p_Te.X[:,0]<0.2, p_Te.y<1)) # Te inside r/a=0.2 > 1keV
    p_ne.remove_points(np.logical_and(p_ne.X[:,0]<0.2, p_ne.y<0.5)) # ne inside r/a=0.2 > 5e19

    # clean up points with points smaller than errorbars
    p_ne.remove_points(p_ne.y < p_ne.err_y)
    p_Te.remove_points(p_Te.y < p_Te.err_y)

    return p_ne, p_Te



def postfit_filter(p_ne, p_Te, ne, Te, roa_kp):

    # impose positivity
    ne[ne<np.nanmin(p_ne.y)] = np.nanmin(p_ne.y)
    Te[Te<np.nanmin(p_Te.y)] = np.nanmin(p_Te.y)

    # eliminate points that are more than 3 sigma away and fit again
    p_ne.remove_points(p_ne.X[:,0]>1.1) # remove artificial point
    chi_ne = (p_ne.y - interp1d(roa_kp, ne)(p_ne.X[:,0]))/p_ne.err_y
    p_ne.remove_points(chi_ne>3)
    p_Te.remove_points(p_Te.X[:,0]>1.1) # remove artificial point
    chi_Te = (p_Te.y - interp1d(roa_kp, Te)(p_Te.X[:,0]))/p_Te.err_y
    p_Te.remove_points(chi_Te>3)

    return p_ne, p_Te
   
 
def create_pe(p_ne, p_Te):

    # use whichever has less as template for pressure
    p_template_X = p_Te.X[:,0] if len(p_Te.X) <= len(p_ne.X) else p_ne.X[:,0]
    p_check_X = p_ne.X[:,0] if len(p_Te.X) <= len(p_ne.X) else p_Te.X[:,0]
       
    p_template_y = p_Te.y if len(p_Te.X) <= len(p_ne.X) else p_ne.y
    p_check_y = p_ne.y if len(p_Te.X) <= len(p_ne.X) else p_Te.y
        
    p_template_err_y = p_Te.err_y if len(p_Te.X) <= len(p_ne.X) else p_ne.err_y
    p_check_err_y = p_ne.err_y if len(p_Te.X) <= len(p_ne.X) else p_Te.err_y
        
    template_inds, check_inds = [], []
    template_avg, check_avg = [], []

    for xx in p_template_X:
        if xx in p_check_X:
            template_list = list(np.where(p_template_X == xx)[0])
            template_inds += template_list
            check_list = list(np.where(p_check_X == xx)[0])
            check_inds += check_list 
        
        num_to_avg = np.abs(len(template_list) - len(check_list))
        if len(template_list) > len(check_list):
            if template_list[:num_to_avg+1] not in template_avg:
                template_avg.append(template_list[:num_to_avg+1])
                check_avg.append(check_list)
        if len(check_list) > len(template_list):
            if check_list[:num_to_avg+1] not in check_avg:
                check_avg.append(check_list[:num_to_avg+1])
                template_avg.append(template_list)

    # remove double counts
    template_inds = list(set(template_inds))
    check_inds = list(set(check_inds))

    unlist_templatelist = [ll for subtemp in check_avg for ll in subtemp]
    unlist_checklist = [ll for subcheck in check_avg for ll in subcheck]
    
    template_X = [p_template_X[tX] for tX in template_inds if tX not in unlist_templatelist]
    for subtemp in template_avg: # for some reason list comprehension doesn't work
        template_X += [np.mean(p_template_X[subtemp])]
    check_X = [p_check_X[tX] for tX in template_inds if tX not in unlist_checklist]
    for subcheck in check_avg:
        check_X += [np.mean(p_check_X[subcheck])]
    
    template_y = [p_template_y[ty] for ty in template_inds if ty not in unlist_templatelist]
    for subtemp in template_avg:
        template_y += [np.mean(p_template_y[subtemp])]
    check_y = [p_check_y[ty] for ty in template_inds if ty not in unlist_checklist]
    for subcheck in check_avg:
        check_y += [np.mean(p_check_y[subcheck])]
    
    template_err_y = [p_template_err_y[tyy] for tyy in template_inds if tyy not in unlist_templatelist]
    for subtemp in template_avg:
        template_err_y += [np.mean(p_template_err_y[subtemp])]
    check_err_y = [p_check_err_y[tyy] for tyy in template_inds if tyy not in unlist_checklist]
    for subcheck in check_avg:
        check_err_y += [np.mean(p_check_err_y[subcheck])]
    
    sort_template = np.argsort(template_X)
    sort_check = np.argsort(check_X)
 
    p_pe_X = np.array(template_X)[sort_template]
    p_pe_y = np.array(template_y)[sort_template]*np.array(check_y)[sort_check]
    p_pe_err_y = np.sqrt(np.array(template_err_y)[sort_template]**2 + np.array(check_err_y)[sort_check]**2)

    from profiletools import BivariatePlasmaProfile as BPP
    p_pe = BPP(X_dim=1,
                X_units=[''],
                y_units='10^{20}$ keV m^{-3}$',
                X_labels=['$r/a$'],
                y_label=r'$p_e$')
    
    p_pe.add_data(p_pe_X, p_pe_y, err_y=p_pe_err_y)

    return p_pe



def fetch_edge_probes(shot, time, Te_sep_eV, geqdsk=None, roa_min=0.995, roa_max=1.05, probes=['A'], max_ne_cm3=5e13, shift_probes=True):
    '''Load data for the ASP and FSP probes on Alcator C-Mod. 
    Time in seconds.

    rhop_min and rhop_max are used to subselect the radial range of the data.

    This function returns profiletools data structures.
    See https://profiletools.readthedocs.io/en/latest/#

    Parameters (TODO: complete description)
    ----------
    shot
    time
    Te_sep_eV
    geqdsk
    roa_min
    roa_max
    probes : list, subselection of ['A','F']
    max_ne_cm3 : maximum density to be permitted. All data points with density greater than this value will be eliminated.
        Default is to have this value extremely high (basically, not active.

    Returns
    -------
    p_ne, p_Te : profiletools object containing combined data from the fetched probes.

    '''
    if geqdsk is None:
        # the geqdsk is used only when using the Aurora radial coordinate transformations
        geqdsk = get_geqdsk_cmod(
            shot, time*1e3, gfiles_loc = '/home/millerma/lya/gfiles/')

    import afsp_probes
    
    has_A = False
    has_F = False
    if 'A' in probes:
      try:
          # if available, add A-Side Probe (ASP) data

          # check raw data first
          out_asp = afsp_probes.get_clean_data(shot, time, geqdsk, probe='A', plot=False)
          roa_asp, roa_unc_asp, t_range, ne_prof_asp, ne_unc_prof_asp, Te_prof_asp, Te_unc_prof_asp, ax = out_asp

          # these variables defined for conditional statements later
          ne_adb = []
          offset_adb = None
          
          # check edge db as well
          out_adb = afsp_probes.get_edgedb_data(shot, t_range, geqdsk, probe='A', plot=False)
          roa_adb, ne_adb, Te_adb, offset_adb = out_adb

          # if edge db populated, take that data
          if any(ne_adb):
              
              ne_not_none = np.where(ne_adb != None)
              Te_not_none = np.where(Te_adb != None)

              if offset_adb is None:
                  offset_adb = np.abs(roa_adb[ne_not_none][0] - roa_asp[0]) # sometimes offset not populated even though roa are different

              print(f'Using probe profiles from edge database')
              print(f'Shifting of ASP data from database, in r/a units: {offset_adb:.4f}')

              # edge db doesn't contain uncertainties so interpolate
              interp_shift = roa_asp[0] - roa_adb[ne_not_none][0]
              ne_unc_prof_adb = interp1d(roa_asp, ne_unc_prof_asp, bounds_error=None, fill_value='extrapolate')(roa_adb+interp_shift) 
              Te_unc_prof_adb = interp1d(roa_asp, Te_unc_prof_asp, bounds_error=None, fill_value='extrapolate')(roa_adb+interp_shift) 
             
              # find where None values are in arrays and set old values to new values
              ne_roa_asp, Te_roa_asp = roa_adb[ne_not_none], roa_adb[Te_not_none] # sometimes they're different...
              ne_prof_asp, Te_prof_asp = ne_adb[ne_not_none], Te_adb[Te_not_none]
              ne_unc_prof_asp, Te_unc_prof_asp = ne_unc_prof_adb[ne_not_none], Te_unc_prof_adb[Te_not_none]

              ne_unc_prof_asp = np.maximum(ne_unc_prof_asp, 0)
              Te_unc_prof_asp = np.maximum(Te_unc_prof_asp, 0)

          else:
              print('No probe data in edge database - using data in tree')

          # probe data is returned in SI units, change to units of 1e20m^-3 and keV
          ne_prof_asp /= 1e20
          ne_unc_prof_asp /= 1e20
          Te_prof_asp /= 1e3
          Te_unc_prof_asp /= 1e3


          # mask out data points outside of given range
          mask_ne_asp = (ne_roa_asp>roa_min)&(ne_roa_asp<roa_max)
          mask_Te_asp = (Te_roa_asp>roa_min)&(Te_roa_asp<roa_max)

          ne_ASP_X = np.ones((len(ne_roa_asp),1))
          Te_ASP_X = np.ones((len(Te_roa_asp),1))
          ne_ASP_X[:,0] = ne_roa_asp
          Te_ASP_X[:,0] = Te_roa_asp

          p_ne_ASP = profiletools.BivariatePlasmaProfile(X_dim=1, X_units='', y_units='$10^{20}$ m$^{-3}$',
                                                         X_labels=r'$r/a$', y_label=r'$n_e$, ASP')
          p_ne_ASP.abscissa = 'r/a'
          p_ne_ASP.shot = shot
          p_ne_ASP.t_min = t_range[0]; p_ne_ASP.t_max = t_range[1]
          channels = range(0, len(ne_prof_asp[mask_ne_asp]))
          p_ne_ASP.add_data(ne_ASP_X[mask_ne_asp], ne_prof_asp[mask_ne_asp], channels={0: channels}, err_y=ne_unc_prof_asp[mask_ne_asp])
          p_ne_ASP.remake_efit_tree() # needed to add this for some reason on 11/1/22

          p_Te_ASP = profiletools.BivariatePlasmaProfile(X_dim=1, X_units='', y_units='keV',
                                                         X_labels=r'$r/a$', y_label=r'$T_e$, ASP')
          p_Te_ASP.abscissa = 'r/a'
          p_Te_ASP.shot = shot
          p_Te_ASP.t_min = t_range[0]; p_Te_ASP.t_max = t_range[1]
          channels = range(0, len(Te_prof_asp[mask_Te_asp]))
          p_Te_ASP.add_data(Te_ASP_X[mask_Te_asp], Te_prof_asp[mask_Te_asp], channels={0: channels}, err_y=Te_unc_prof_asp[mask_Te_asp])
          p_Te_ASP.remake_efit_tree() # needed to add this for some reason on 11/1/22
          
          has_A = True
      except Exception:
          print('ASP fetch failed')
          pass

    ######
    if 'F' in probes:
      try:
          # if available, add F-Side Probe (FSP) data

          # check raw data first
          out = afsp_probes.get_clean_data(shot, time, geqdsk, probe='F', plot=False)
          roa_fsp, roa_fsp_unc, t_range, ne_prof_fsp, ne_unc_prof_fsp, Te_prof_fsp, Te_unc_prof_fsp, ax = out

          # define variables for conditionals
          ne_fdb = []
          offset_fdb = None

          # check edge db as well
          out_fdb = afsp_probes.get_edgedb_data(shot, t_range, geqdsk, probe='F', plot=False)
          roa_fdb, ne_fdb, Te_fdb, offset_fdb = out_fdb
          

          # if edge db populated, take that data
          if any(ne_fdb):

              ne_not_none = np.where(ne_fdb != None)
              Te_not_none = np.where(Te_fdb != None)
              
              if offset_fdb is None:
                  offset_fdb = np.abs(roa_fdb[ne_not_none][0] - roa_fsp[0]) # sometimes offset not populated even though roa are different
              
              print(f'Using probe profiles from edge database')
              print(f'Shifting of FSP data from database, in r/a units: {offset_fdb:.4f}')

              # edge db doesn't contain uncertainties so interpolate
              interp_shift = roa_fsp[0] - roa_fdb[ne_not_none][0]
              ne_unc_prof_fdb = interp1d(roa_fsp, ne_unc_prof_fsp, bounds_error=None, fill_value='extrapolate')(roa_fdb+interp_shift) 
              Te_unc_prof_fdb = interp1d(roa_fsp, Te_unc_prof_fsp, bounds_error=None, fill_value='extrapolate')(roa_fdb+interp_shift) 
             
              # find where None values are in arrays and set old values to new values
              ne_roa_fsp, Te_roa_fsp = roa_fdb[ne_not_none], roa_fdb[Te_not_none] # sometimes they're different...
              ne_prof_fsp, Te_prof_fsp = ne_fdb[ne_not_none], Te_fdb[Te_not_none]
              ne_unc_prof_fsp, Te_unc_prof_fsp = ne_unc_prof_fdb[ne_not_none], Te_unc_prof_fdb[Te_not_none]

              ne_unc_prof_fsp = np.maximum(ne_unc_prof_fsp, 0)
              Te_unc_prof_fsp = np.maximum(Te_unc_prof_fsp, 0)
              
          else:
              print('No probe data in edge database - using data in tree')


          # probe data is returned in SI units, change to units of 1e20m^-3 and keV
          ne_prof_fsp /= 1e20
          ne_unc_prof_fsp /= 1e20
          Te_prof_fsp /= 1e3
          Te_unc_prof_fsp /= 1e3

          # mask out data points outside of given range
          mask_ne_fsp = (ne_roa_fsp>roa_min)&(ne_roa_fsp<roa_max)
          mask_Te_fsp = (Te_roa_fsp>roa_min)&(Te_roa_fsp<roa_max)

          ne_FSP_X = np.ones((len(ne_roa_fsp),1))
          ne_FSP_X[:,0] = ne_roa_fsp
          Te_FSP_X = np.ones((len(Te_roa_fsp),1))
          Te_FSP_X[:,0] = Te_roa_fsp

          p_ne_FSP = profiletools.BivariatePlasmaProfile(X_dim=1, X_units='', y_units='$10^{20}$ m$^{-3}$',
                                                         X_labels=r'$r/a$', y_label=r'$n_e$, FSP')
          p_ne_FSP.abscissa = 'r/a'
          p_ne_FSP.shot = shot
          p_ne_FSP.t_min = t_range[0]; p_ne_FSP.t_max = t_range[1]        
          channels = range(0, len(ne_prof_fsp[mask_ne_fsp]))
          p_ne_FSP.add_data(ne_FSP_X[mask_ne_fsp], ne_prof_fsp[mask_ne_fsp], channels={0: channels}, err_y=ne_unc_prof_fsp[mask_ne_fsp])
          p_ne_FSP.remake_efit_tree() # needed to add this for some reason on 11/1/22

          p_Te_FSP = profiletools.BivariatePlasmaProfile(X_dim=1, X_units='', y_units='keV',
                                                         X_labels=r'$\sqrt{\psi_n}$', y_label=r'$T_e$, FSP')
          p_Te_FSP.abscissa = 'r/a'
          p_Te_FSP.shot = shot
          p_Te_FSP.t_min = t_range[0]; p_Te_FSP.t_max = t_range[1]
          channels = range(0, len(Te_prof_fsp[mask_Te_fsp]))
          p_Te_FSP.add_data(Te_FSP_X[mask_Te_fsp], Te_prof_fsp[mask_Te_fsp], channels={0: channels}, err_y=Te_unc_prof_fsp[mask_Te_fsp])
          p_Te_FSP.remake_efit_tree() # needed to add this for some reason on 11/1/22
          
          has_F = True

      except Exception:
          print('FSP fetch failed')
          pass
        

    def probe_func(x,a,k,b):
        ''' Exponential function used to fit probe data'''
        return a*np.exp(-k*(x-1))+b


    # extended range for fitting -- in practice though, shifting should be <0.005 units typically
    x_out = np.linspace(0.96, 1.05, 300)

    # if offset_[a/f]db is not None, probes have already been shifted - shift_probes flag might be turned off if we try to match pressure instead
    if has_A and shift_probes and offset_adb is None:
        
        # fit exponential over extended range near the LCFS for Te
        popt_Te, pcov_Te = curve_fit(probe_func,p_Te_ASP.X[:,0], p_Te_ASP.y,
                                     p0=[0.06,1e2,0.03]) # important to provide a decent guess
        Te_keV_ASP_fit = probe_func(x_out,*popt_Te)

        try:
            # now, based on this exponential, find location of the LCFS        
            xSep_ASP = interp1d(Te_keV_ASP_fit, x_out, bounds_error=True)(Te_sep_eV*1e-3)
            print(f'Shifting of ASP data from fit, in r/a units: {1-xSep_ASP:.4f}')

            # shift both ne and Te data
            p_Te_ASP.X += 1 - xSep_ASP
            p_ne_ASP.X += 1 - xSep_ASP
            
        except ValueError:
            print('Shifting of ASP data failed, likely probe did not go very far into the plasma')

            has_A=False


    if has_F and shift_probes and offset_fdb is None:
        # fit exponential over extended range near the LCFS for Te
        popt_Te, pcov_Te = curve_fit(probe_func,p_Te_FSP.X[:,0], p_Te_FSP.y,
                                     p0=[0.06,1e2,0.03]) # important to provide a decent guess
        Te_keV_FSP_fit = probe_func(x_out,*popt_Te)

        try:
            # now, based on this exponential, find location of the LCFS        
            xSep_FSP = interp1d(Te_keV_FSP_fit, x_out, bounds_error=True)(Te_sep_eV*1e-3)
            print(f'Shifting of FSP data from fit, in r/a units: {1-xSep_FSP:.4f}')
            
            # shift both ne and Te data
            p_Te_FSP.X += 1 - xSep_FSP
            p_ne_FSP.X += 1 - xSep_FSP

        except ValueError:
            print('Shifting of FSP data failed, likely probe did not go very far into the plasma')

            has_F=False

    import copy
    
    # combine data from ASP and FSP probes if both are available
    if has_A and has_F:
        p_ne = copy.deepcopy(p_ne_ASP)
        p_ne.add_profile(p_ne_FSP)
        p_Te = copy.deepcopy(p_Te_ASP)
        p_Te.add_profile(p_Te_FSP)
    elif has_A and not has_F:
        p_ne = copy.deepcopy(p_ne_ASP)
        p_Te = copy.deepcopy(p_Te_ASP)
    elif has_F and not has_A:
        p_ne = copy.deepcopy(p_ne_FSP)
        p_Te = copy.deepcopy(p_Te_FSP)
    else:
        p_ne = None
        p_Te = None

    # multiply probe ne by some factor to match TS
    #mult_factor = 2

    #p_ne.y*=mult_factor
    #p_ne.y*=mult_factor

    # remove probe points for r/a < 0.99 (typically unreliable)
    if p_ne is not None:
        p_ne.remove_points(p_ne.X[:,0] < 0.99)
        p_Te.remove_points(p_Te.X[:,0] < 0.99)

        # remove data corresponding to likely emissive probe behavior
        p_Te.remove_points(p_ne.y*1e14 > max_ne_cm3)
        p_ne.remove_points(p_ne.y*1e14 > max_ne_cm3)

        if p_ne.X is None:
            raise ValueError('All points from probe ne were removed! max_ne_cm3 is likely too small.')
  
    return p_ne, p_Te


def filter_TS(ne_TS, Te_TS, ne_SP, Te_SP, cutoff):
    
    if ne_SP is not None:

        # make sure X array is increasing
        order_SP = np.argsort(Te_SP.X[:,0])

        filt = np.where(Te_SP.y[order_SP] < cutoff)[0]
        x_filt = Te_SP.X[order_SP][filt[0],0]
        ne_TS.remove_points(ne_TS.X[:,0] > x_filt)
        Te_TS.remove_points(Te_TS.X[:,0] > x_filt)

    return ne_TS, Te_TS, ne_SP, Te_SP


def filter_SP(ne_TS, Te_TS, ne_SP, Te_SP, cutoff):

    # make sure X array is increasing
    order = np.argsort(Te_TS.X[:,0])

    filt = np.where(Te_TS.y[order] > cutoff)[0]
    x_filt = Te_TS.X[order][filt[-1],0]
    ne_SP.remove_points(ne_SP.X[:,0] < x_filt)
    Te_SP.remove_points(Te_SP.X[:,0] < x_filt)

    return ne_TS, Te_TS, ne_SP, Te_SP


def match_pressure(ne_TS, Te_TS, ne_SP, Te_SP, plot=False):

    ne_SP_sep = interp1d(ne_SP.X[:,0], ne_SP.y)(1)
    Te_SP_sep = interp1d(Te_SP.X[:,0], Te_SP.y)(1)
    # calculate electron pressure at roa = 1 as measured by SP
    pe_SP_sep = ne_SP_sep*Te_SP_sep

    # find where the SP pressure matches the TS pressure
    # need to interpolate Te onto ne
    Te_TS_intp = interp1d(Te_TS.X[:,0], Te_TS.y, fill_value='extrapolate')(ne_TS.X[:,0])
    pe_TS = ne_TS.y*Te_TS_intp
    roa_TS = ne_TS.X[:,0]

    roa_TS_match = interp1d(pe_TS, roa_TS, fill_value='extrapolate')(pe_SP_sep)

    # shift SP profiles
    SP_shift = roa_TS_match - 1
    ne_SP.X[:,0] += SP_shift
    Te_SP.X[:,0] += SP_shift
            
    print(f'Shifting of ASP data, in r/a units: {SP_shift:.4f}')

    if plot:

        # calculate overall pressure SP profile to plot in same way as TS

        Te_SP_intp = interp1d(Te_SP.X[:,0], Te_SP.y, fill_value='extrapolate')(ne_SP.X[:,0])
        pe_SP = ne_SP.y*Te_SP_intp
        roa_SP = ne_SP.X[:,0] - SP_shift

        fig, ax = plt.subplots()
        ax.plot(roa_TS, pe_TS, 'o')
        ax.plot(roa_SP, pe_SP, 'o')
        ax.plot(ne_SP.X[:,0], pe_SP, 'o')
        plt.xlim([0.86,1.16])
        plt.ylim([-0.01, 0.27])
        plt.legend(['TS', 'SP before shift', 'SP after shift'])

    return ne_SP, Te_SP


# def fav_vs_unfav(shot,time,geqdsk = None):
#     '''Determine whether grad-B field direction is favorable or unfavorable.
    
#     This function ignores the possibility of having a double null, use with care!
#     The active x-point is taken to be the point along the LCFS that is furthest from Z=0.

#     '''
#     if geqdsk is None:
#         geqdsk = get_geqdsk_cmod(shot,time*1e3)
#         geqdsk = get_geqdsk_cmod(shot,time*1e3)  # repeat to make sure it's loaded...
    
#     # find ion grad(B)-drift direction (determined by B field dir, since radial grad(B) is always inwards )
#     magTree = MDSplus.Tree('magnetics',shot)
#     nodeBt = magTree.getNode('\magnetics::Bt')
#     Bt = nodeBt.data()
#     time_Bt = nodeBt.dim_of().data()
#     tidx = np.argmin(np.abs(time_Bt - time)) 
#     gradB_drift_up = False if Bt[tidx]<0 else True
    
#     # find whether shot is USN or LSN -- assume not DN...
#     maxZ = np.max(geqdsk['ZBBBS'])
#     minZ = np.min(geqdsk['ZBBBS'])
    
#     #  pretty sure that the X-point is where the LCFS is furthest from the magnetic axis
#     USN = True if np.abs(maxZ)==np.max([np.abs(maxZ),np.abs(minZ)]) else False
    
#     # favorable or unfavorable grad-B drift direction?
#     favorable = (gradB_drift_up and USN) or (gradB_drift_up==False and USN==False)
    
#     return gradB_drift_up, USN, favorable

        
def get_vol_avg_pressure(shot,time,rhop,ne,Te,geqdsk=None):
    ''' Calculate volume-averaged pressure given some ne,Te radial profiles.

    ne must be in cm^-3 units and Te in eV.
    '''
    # find volume-averaged pressure    
    p_Pa = (ne*1e6) * (Te*q_electron)
    p_atm = p_Pa/101325.  # conversion factor between Pa and atm

    # load geqdsk dictionary
    if geqdsk is None:
        geqdsk = get_geqdsk_cmod(shot,time*1e3)
    
    # find volume average within LCFS
    indLCFS = np.argmin(np.abs(rhop-1.0))
    p_Pa_vol_avg = aurora.vol_average(p_Pa[:indLCFS], rhop[:indLCFS], method='omfit',geqdsk = geqdsk)[-1]
    #p_atm_vol_avg = p_Pa_vol_avg/101325.

    return p_Pa_vol_avg



def get_geqdsk_cmod(shot, time_ms, gfiles_loc = '/home/sciortino/EFIT/gfiles/'):
    ''' Get a geqdsk file in omfit_eqdsk format by loading it from disk, if available, 
    or from MDS+ otherwise.  

    This function tries to first load a EFIT20, if available.

    time must be in ms!

    Currently, the omfit_eqdsk class struggles to connect to MDS+ sometimes. To work 
    around this problem, this function uses a loop and try-except statements to try to 
    load multiple times until succeeding.
    '''
    time_ms = np.floor(time_ms)   # TODO: better to do this outside!!
    file_name=f'g{shot}.{str(int(time_ms)).zfill(5)}'

    def fetch_and_move():
        try:
            # try to fetch EFIT20 first
            geqdsk = omfit_eqdsk.OMFITgeqdsk('').from_mdsplus(
                device='CMOD',shot=shot, time=time_ms, SNAPfile='EFIT20',
                fail_if_out_of_range=True,time_diff_warning_threshold=20
            )
        except:
            # if EFIT20 is not available, look for default ANALYSIS EFIT
            geqdsk = omfit_eqdsk.OMFITgeqdsk('').from_mdsplus(
                device='CMOD',shot=shot, time=time_ms, SNAPfile='ANALYSIS',
                fail_if_out_of_range=True,time_diff_warning_threshold=20
            )
    
        geqdsk.save(raw=True)
        shutil.move(file_name, gfiles_loc+file_name)


    def attempt_loading():
        if os.path.exists(gfiles_loc + file_name):
            # fetch local g-file if available
            try:
                geqdsk = omfit_eqdsk.OMFITgeqdsk(gfiles_loc + file_name)
                kk = geqdsk.keys()  # quick test
            except:
                geqdsk = fetch_and_move()
        else:
            geqdsk = fetch_and_move()

        return geqdsk
    
    # sometimes there are issues loading and this must be enforced as follows
    # this may become redundant in the future when this bug is solved...
    for ijk in np.arange(10): # 10 attempts
        try:
            geqdsk = attempt_loading()
            geqdsk.load()  # for safety
            geqdsk['fluxSurfaces'].load()  # for safety
            break
        except Exception as e:
            pass
    
    return geqdsk

    
def get_Greenwald_frac(shot, tmin,tmax, roa, ne, Ip_MA, a_m=0.22, geqdsk=None):
    ''' Calculate Greenwald density fraction by normalizing volume-averaged density.

    INPUTS
    ------
    shot : int, shot number
    tmin and tmax: floats, time window (in [s]) to fetch equilibrium.
    ne: 1D array-like, expected as time-independent. Units of 1e20 m^-3.
    Ip_MA: float, plasma current in MA. 
    a_m : minor radius in units of [m]. Default of 0.69 is for C-Mod. 

    OUTPUTS
    -------
    n_by_nG : float
        Greenwald density fraction, defined with volume-averaged density.
    '''
    if geqdsk is None:
        tmin *= 1000.  # change to ms
        tmax *= 1000.  # change to ms
        time = (tmax+tmin)/2.
        geqdsk = get_geqdsk_cmod(shot,time)
    
    # find volume average within LCFS
    rhop = aurora.rad_coord_transform(roa,'r/a','rhop', geqdsk)

    indLCFS = np.argmin(np.abs(rhop-1.0))
    n_volavg = aurora.vol_average(ne[:indLCFS], rhop[:indLCFS], geqdsk=geqdsk)[-1]

    # Greenwald density
    n_Gw = Ip_MA/(np.pi * a_m**2)   # units of 1e20 m^-3, same as densities above

    # Greenwald fraction:
    f_gw = n_volavg/n_Gw

    return f_gw



def get_CMOD_gas_fueling(shot, plot=False):
    '''Load injected gas amounts and give a grand total in Torr-l.
    Translated from gas_input2_ninja.dat scope. 
    '''

    _c_side = smooth(omfit_mds.OMFITmdsValue(server='CMOD',shot=shot,treename='cmod',
                                             TDI='\\plen_cside').data()[0,:],31)
    _t = omfit_mds.OMFITmdsValue(server='CMOD',shot=shot,treename='cmod',
                                 TDI='dim_of(\\plen_cside)').data()
    _b_sideu = smooth(omfit_mds.OMFITmdsValue(server='CMOD',shot=shot,treename='cmod',
                                              TDI='\\plen_bsideu').data()[0,:],31)
    _b_top = smooth(omfit_mds.OMFITmdsValue(server='CMOD',shot=shot,treename='cmod',
                                            TDI='\\plen_btop').data()[0,:],31)

    ninja = True
    try:
        plen_bot_time = omfit_mds.OMFITmdsValue(server='CMOD',shot=shot,treename='edge',
                                            TDI='\edge::gas_ninja.plen_bot').dim_of(0)
        plen_bot = smooth(omfit_mds.OMFITmdsValue(server='CMOD',shot=shot,treename='edge',
                                              TDI='\edge::gas_ninja.plen_bot').data()[0,:],31)
    except:
        ninja = False
        print('No gas fueling frlom ninja system')

    # only work with quantities within [0,2]s interval
    ind0 = np.argmin(np.abs(_t)); ind1 = np.argmin(np.abs(_t-2.0))

    time = _t[ind0:ind1]
    c_side = _c_side[ind0:ind1]
    b_sideu = _b_sideu[ind0:ind1]
    b_top = _b_top[ind0:ind1]

    # ninja system is on a different time base than the other measurements
    ninja2 = interp1d(plen_bot_time, plen_bot, bounds_error=False)(time) if ninja else 0
    
    gas_tot = c_side + b_sideu + b_top + ninja2

    if plot:
        fig,ax = plt.subplots()
        ax.plot(time, gas_tot, label='total')
        ax.plot(time, c_side, label='c-side')
        ax.plot(time, b_sideu, label='b-side u')
        ax.plot(time, b_top, label='b-top')
        ax.plot(time, ninja2, label='ninja2')
        ax.legend(loc='best').set_draggable(True)
        ax.set_xlabel('time [s]')
        ax.set_ylabel('Total injected gas [Torr-l]')
        
    return time, gas_tot

    
def get_Lya_data(shot=1080416024, systems=['LYMID'], plot=True):
    ''' Get Ly-alpha data for C-Mod from any (or all) of the systems:
    ['LYMID','WB1LY','WB4LY','LLY','BPLY']
    '''

    bdata = {} # BRIGHT node
    edata = {} # EMISS node

    if systems=='all':
        systems=['LYMID','WB1LY','WB4LY','LLY','BPLY']
        
    if plot:
        fig, ax = plt.subplots(1,2, figsize=(13,8))
        ls = ['-','--','-.',':','--']

    for ss,system in enumerate(systems):
        fetched_0=True; fetched_1=True
        try:
            bdata[system] = fetch_bright(shot, system)
            
            if plot:
                for ch in np.arange(bdata[system].values.shape[1]):
                    ax[0].plot(bdata[system].time, bdata[system].values[:,ch],
                               label=system+', '+str(ch), ls=ls[ss])
        except:
            print('Could not fetch C-Mod Ly-alpha BRIGHT data from system '+system)
            fetched_0=False
            pass
        
        try:
            edata[system] = fetch_emiss(shot, system)
            if plot:
                for ch in np.arange(edata[system].values.shape[1]):
                    ax[1].plot(edata[system].time, edata[system].values[:,ch],
                               label=system+', '+str(ch), ls=ls[ss])
        except:
            print('Could not fetch C-Mod Ly-alpha EMISS data from system '+system)
            fetched_1=False
            pass
        
        if plot:
            ax[0].set_xlabel('time [s]'); ax[1].set_xlabel('time [s]')
            if fetched_0: ax[0].set_ylabel(r'Brightness [$'+str(bdata[system].units)+'$]')
            if fetched_1: ax[1].set_ylabel(r'Emissivity [$'+str(edata[system].units)+'$]')
            ax[0].legend(); ax[1].legend()
            
    return bdata,edata


def calc_nu_star(ne, Te):

    ne = ne*1e20
    Te = Te*1e3

    m_i = 1.67e-27
    m_e = 9.11e-31
    mu = m_i/m_e
    q_e = 1.602e-19
    e0 = 8.85e-12
    a = 0.22
    R0 = 0.68
    safety_factor = 2 ### how to get safety factor??
    
    v_th = np.sqrt(Te*q_e/m_e)
    nu = 1/(4*np.pi)*q_e**4*ne*15/(e0**2*m_e**(1/2)*(Te*q_e)**(3/2))
    epsilon = a/R0
    nu_star = nu/epsilon**(3/2)/(v_th/safety_factor*R0)

    return nu_star


def identify_mode(shot):

    lrf07_array = np.genfromtxt('lyman_lmode_fy07_test.txt',skip_header=2)
    loh07_array = np.genfromtxt('lyman_ohmic_fy07_test.txt',skip_header=2)
    h07_array = np.genfromtxt('lyman_hmode_fy07_test.txt',skip_header=2)
    
    lrf08_array = np.genfromtxt('lyman_lmode_fy08_test.txt',skip_header=2)
    loh08_array = np.genfromtxt('lyman_ohmic_fy08_test.txt',skip_header=2)
    h08_array = np.genfromtxt('lyman_hmode_fy08_test.txt',skip_header=2)
    i08_array = np.genfromtxt('lyman_imode_fy08_test.txt',skip_header=2)

    l_modes = [int(s) for s in  lrf07_array[:,1]]
    l_modes += [int(s) for s in  loh07_array[:,1]]
    l_modes += [int(s) for s in  lrf08_array[:,1]]
    l_modes += [int(s) for s in  loh08_array[:,1]]
    h_modes = [int(s) for s in  h07_array[:,1]]
    h_modes += [int(s) for s in  h08_array[:,1]]
    i_modes = [int(s) for s in  i08_array[:,1]]

    if shot in l_modes:
        mode = 'L'
    elif shot in h_modes:
        mode = 'H'
    elif shot in i_modes:
        mode = 'I'
    else:
        print('Shot not in database')
        mode = input('Specify mode if known (H/L/I), else L-mode will be assumed:')
        if mode == 'N': mode = 'L' 
        return mode

    return mode


def get_fit_gradient(y, c, x, fit_type, grad_type='analytic', reg=None, plot=True):

    if grad_type == 'analytic':

        if fit_type == 'osborne': grad = mtanh_fitting.Osbourne_Tanh_gradient(x, c, reg=reg)
        if fit_type == 'superfit': grad = mtanh_fitting.mtanh_profile_gradient(x, edge=c[0], ped=c[1], core=c[2], expin=c[3], expout=c[4], widthp=c[5], xphalf=c[6])

    else:
        grad = np.gradient(y)

    if plot:

        fig, ax = plt.subplots(2, sharex=True)
        ax[0].plot(x, y, 'k-')    
        ax[1].plot(x, grad, 'k-')

    return grad


class kinetic_profile:

    # just a little class to store results from fits so that we can do monte carlo error estimation more easily
    
    def __init__(self, x, y, grad_y):

        self.x = x
        self.y = y
        self.grad_y = grad_y

#def perturb_mc_err(num_mc, p_ne, p_Te, p_pe, ne, Te, pe, roa_kp, maxfev, ped_width, reg=[], remove_outliers=False, verbose=True, plot=False):
    def perturb_mc_err(self, raw, num_mc, maxfev, ped_width, reg, ne=True, verbose=True, plot=False):
        
        err = raw.err_y

        if plot:
    
            fig, ax = plt.subplots()
            ax.errorbar(raw.X[:-1,0], raw.y[:-1], raw.err_y[:-1], fmt='.')
            ax.set_ylabel('$n_{e}$ (10$^{20}$m$^{-3})$')
            ax.set_xlabel('r/a')


        max_iter = np.max(num_mc) # get biggest one
    
        if type(num_mc) == int: num_mc = np.array([num_mc]) # if pass in int, convert to array
            
        # these arrays hold standard deviations and mean for different number of iterations
        npts_fit = len(self.y) # number of points in fit

        # these arrays hold the overall distributions of ne, Te, and their chisqr
        fit_dist = np.zeros([max_iter, npts_fit])
        grad_fit_dist = np.zeros([max_iter, npts_fit])
        xsqr_dist = np.zeros([max_iter])

        # iterate through MC perturbation of largest number

        for ii in range(max_iter):
            new_ne_pts = raw.y + np.random.normal(scale=raw.err_y)
    
            _out_ne = mtanh_fitting.best_osbourne(raw.X[:,0], new_ne_pts, vals_unc=raw.err_y, x_out=self.x, maxfev=maxfev, reg=reg, ped_width=ped_width, ne=ne)

            fit_dist[ii], popt, perr, xsqr_dist[ii] = _out_ne

            if verbose: 
                print('MC iteration {}'.format(ii), 'XSQR = {:.3f}'.format(xsqr_dist[ii]))

            grad_fit_dist[ii] = get_fit_gradient(fit_dist[ii], popt, self.x, 'osborne', grad_type='analytic', reg=reg, plot=False)
    
    
            if plot: ax.plot(self.x, fit_dist[ii], 'k-')


        # now sort chisqr distributions and get rid of top 10%
        num_mc_10p = np.int(np.ceil(num_mc/10))

        xsqr_inds = np.argsort(xsqr_dist)[-num_mc_10p:]
        keep_mask = np.ones(len(xsqr_dist), dtype=bool)
        keep_mask[xsqr_inds] = False

        
        fit_mean = np.zeros([len(num_mc),npts_fit])
        fit_std = np.zeros([len(num_mc),npts_fit])
        grad_fit_mean = np.zeros([len(num_mc),npts_fit])
        grad_fit_std = np.zeros([len(num_mc),npts_fit])

        # now iterate through number of iterations and take standard deviations of all of them
        for ss in range(len(num_mc)):
        
            fit_std[ss] = np.std(fit_dist[:num_mc[ss]][keep_mask], axis=0)
            fit_mean[ss] = np.mean(fit_dist[:num_mc[ss]][keep_mask], axis=0)
            grad_fit_std[ss] = np.std(grad_fit_dist[:num_mc[ss]][keep_mask], axis=0)
            grad_fit_mean[ss] = np.mean(grad_fit_dist[:num_mc[ss]][keep_mask], axis=0)

            # z_score = 3
            # if remove_outliers:
        
            #     # throw out points in fits that are more than 3std away
            #     ne_norm = np.abs(ne_fit_dist[:num_mc[ss]] - ne_fit_mean[ss])
            #     grad_ne_norm = np.abs(grad_ne_fit_dist[:num_mc[ss]] - grad_ne_fit_mean[ss])
            #     Te_norm = np.abs(Te_fit_dist[:num_mc[ss]] - Te_fit_mean[ss])
            #     grad_Te_norm = np.abs(grad_Te_fit_dist[:num_mc[ss]] - grad_Te_fit_mean[ss])

            #     mask_3std_ne = np.logical_and(ne_norm <= z_score*ne_fit_std[ss], grad_ne_norm <= z_score*grad_ne_fit_std[ss])
            #     ne_fit_std[ss] = np.array([np.std(ne_fit_dist[:,rr][mask_3std_ne[:,rr]]) for rr in range(len(roa_kp))])
            #     grad_ne_fit_std[ss] = np.array([np.std(grad_ne_fit_dist[:,rr][mask_3std_ne[:,rr]]) for rr in range(len(roa_kp))])
        
            #     mask_3std_Te = np.logical_and(Te_norm <= z_score*Te_fit_std[ss], grad_Te_norm <= z_score*grad_Te_fit_std[ss])
            #     Te_fit_std[ss] = np.array([np.std(Te_fit_dist[:,rr][mask_3std_Te[:,rr]]) for rr in range(len(roa_kp))])
            #     grad_Te_fit_std[ss] = np.array([np.std(grad_Te_fit_dist[:,rr][mask_3std_Te[:,rr]]) for rr in range(len(roa_kp))])

        z_score = 2

        conf_up = fit_mean + z_score*fit_std 
        conf_down = fit_mean - z_score*fit_std 
    
        num_out = 0
        in_mask = np.ones(fit_dist.shape[0], dtype=bool)
        outlier_inds = []

        ind_cut = 0 # index from where to begin the confidence interval (0 is full profile, 150 is outer half)

        for sss in range(len(num_mc)):
                for mc in range(max_iter):
   
                    prof_out = np.logical_or(fit_dist[mc][ind_cut:] > conf_up[0,ind_cut:], fit_dist[mc][ind_cut:] < conf_down[0,ind_cut:])

                    if np.sum(prof_out) > 15:
                        num_out += 1
                        outlier_inds.append(mc)
                
        keep_mask[outlier_inds] = False
        removed_inds = np.union1d(outlier_inds, xsqr_inds)

        print('removed inds: {} ({}, {})'.format(len(removed_inds), len(outlier_inds), len(xsqr_inds)))

        if plot:
            conf_up = fit_mean + fit_std 
            conf_down = fit_mean - fit_std 
       
            maxval = np.max(raw.y[raw.X[:,0] > 0.85])
            headroom = 1.05      

            ax.plot(self.x, conf_up[0], 'r-')
            ax.plot(self.x, conf_down[0], 'r-')
        
            ax.set_xlim([0.85,1.1])
            ax.set_ylim([0, maxval*headroom])
        
            # plot chi_square distribution
            fig, ax = plt.subplots()
            ax[0].hist(xsqr_dist, bins=50)
    
    
            # plot and recompute statistics w/o outliers

            fig, ax = plt.subplots()
            ax.errorbar(raw.X[:-1,0], raw.y[:-1], raw.err_y[:-1], fmt='.')
            ax.set_ylabel('$n_{e}$ (10$^{20}$m$^{-3})$')
            ax.set_xlabel('r/a')
        
            for ii in range(max_iter):
                if ii not in ne_outlier_inds:
                    ax.plot(roa_kp, fit_dist[ii], 'k-')
                    ax.plot(roa_kp, conf_up[0], 'r-')
                    ax.plot(roa_kp, conf_down[0], 'r-')
                
            ax.set_xlim([0.85,1.1])
            ax.set_ylim([0, maxval*headroom])
        

        # final fits
        fit_dist_filtered = fit_dist[keep_mask]


        # recompute statistics
        for ss in range(len(num_mc)):
        
            fit_mean[ss] = np.mean(fit_dist[:num_mc[ss]][keep_mask], axis=0)
            fit_std[ss] = np.std(fit_dist[:num_mc[ss]][keep_mask], axis=0)
            grad_fit_mean[ss] = np.mean(grad_fit_dist[:num_mc[ss]][keep_mask], axis=0)
            grad_fit_std[ss] = np.std(grad_fit_dist[:num_mc[ss]][keep_mask], axis=0)

        
        self.fit_dist_prefilt = fit_dist
        self.fit_mean = fit_mean[0] # need to index like this because of earlier test done to find what best MC number is
        self.fit_std = fit_std[0]
        self.grad_fit_mean = grad_fit_mean[0]
        self.grad_fit_std = grad_fit_std[0]
        self.fit_dist = fit_dist_filtered

        return None
   
 
    def save_coefs(self, fit_popt, mean_popt):
    
        self.popt = np.vstack((fit_popt, mean_popt))
    


def fetch_bright(shot,system):

    _bdata = {}
    node = omfit_mds.OMFITmdsValue(server='CMOD', shot=shot, treename='SPECTROSCOPY',
                         TDI='\\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.'+\
                             '{:s}:BRIGHT'.format(system))
    _bdata = xarray.DataArray(
        node.data(), coords={'time':node.dim_of(1),'R':node.dim_of(0)},
        dims=['time','R'],
        attrs={'units': node.units()})

    return _bdata


def fetch_emiss(shot,system):

    _edata = {}
    node = omfit_mds.OMFITmdsValue(server='CMOD', shot=shot, treename='SPECTROSCOPY',
                         TDI='\\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.'+\
                             '{:s}:EMISS'.format(system))
    _edata = xarray.DataArray(
        node.data(), coords={'time':node.dim_of(1),'R':node.dim_of(0)},
        dims=['time','R'],
        attrs={'units': node.units()})

    #print('Emissivity units: ' , node.units())

    return _edata


def fetch_tomo_emiss(shot,system,r_end=0.93,sys_err=5,shift=0):

    _out = tomo.tomoCMOD(shot,system,r_end=r_end,sys_err=sys_err)
    tvec,R_grid,y,y_err,backprojection = _out

    _edata = xarray.DataArray(
        y, coords={'time':tvec,'R':R_grid+shift},
        dims=['time','R',],
        attrs={'units': '$W/m^{3}$'})

    return _edata, y_err


def get_CMOD_1D_geom(shot,time):

    # right gap
    tmp = omfit_mds.OMFITmdsValue(server='CMOD', treename='ANALYSIS', shot=shot,
                        TDI='\\ANALYSIS::TOP.EFIT.RESULTS.A_EQDSK.ORIGHT')
    time_vec = tmp.dim_of(0)
    _gap_R = tmp.data()
    gap_R = _gap_R[time_vec.searchsorted(time)-1]

    # R location of LFS LCFS
    tmp = omfit_mds.OMFITmdsValue(server='CMOD', treename='ANALYSIS', shot=shot,
                        TDI='\\ANALYSIS::TOP.EFIT.RESULTS.G_EQDSK.RBBBS')
    time_vec = tmp.dim_of(0)
    _rbbbs = tmp.data()*1e2 # m --> cm
    rbbbs = _rbbbs[:,time_vec.searchsorted(time)-1]

    Rsep = np.max(rbbbs)

    return Rsep,gap_R




def plot_emiss(edata, shot, time, ax=None):
    ''' Plot emissivity profile '''

    # get Rsep and gap
    Rsep, gap = get_CMOD_1D_geom(shot,time)
    Rwall = Rsep+gap
    print('Rwall,Rsep,gap:',Rwall, Rsep,gap)

    if ax is None:
        fig,ax = plt.subplots()

    tidx = np.argmin(np.abs(edata.time.values - time))
    ax.plot(edata.R.values, edata.values[tidx,:], '.-')  #*100 - Rwall

    ax.set_ylabel(r'emissivity [${:}$]'.format(edata.units))
    ax.set_xlabel(r'R [cm]')
    return ax


def plot_bright(bdata, shot, time,ax=None):
    ''' Plot brightness over chords profile '''

    # get Rsep and gap
    Rsep, gap = get_CMOD_1D_geom(shot,time)
    Rwall = Rsep+gap
    print('Rwall,Rsep,gap:',Rwall, Rsep,gap)

    if ax is None:
        fig,ax = plt.subplots()

    tidx = bdata.time.values.searchsorted(time)-1
    mask = np.nonzero(bdata.values[tidx,:])[0]
    ax.plot(bdata.R.values[mask], bdata.values[tidx,mask], '.-')  #*100-Rwall

    ax.set_ylabel(r'brightness [${:}$]'.format(bdata.units))
    ax.set_xlabel(r'R [cm]')

    return ax

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def get_P_ohmic(shot):
    ''' Get Ohmic power

    Translated/adapted from scopes:
    _vsurf =  deriv(smooth1d(\ANALYSIS::EFIT_SSIBRY,2))*$2pi ;
    _ip=abs(\ANALYSIS::EFIT_AEQDSK:CPASMA);
    _li = \analysis::efit_aeqdsk:ali;
    _L = _li*6.28*67.*1.e-9;
    _vi = _L*deriv(smooth1d(_ip,2));
    _poh=_ip*(_vsurf-_vi)/1.e6
    '''

    # psi at the edge:
    ssibry_node = OMFITmdsValue(server='CMOD', shot=shot, treename='ANALYSIS',TDI='\\analysis::efit_ssibry')
    time = ssibry_node.dim_of(0)
    ssibry = ssibry_node.data()
    
    # total voltage associated with magnetic flux inside LCFS
    vsurf = np.gradient(smooth(ssibry,5),time) * 2 * np.pi

    # calculated plasma current
    ip_node= OMFITmdsValue(server='CMOD', shot=shot, treename='ANALYSIS',TDI='\\analysis::EFIT_AEQDSK:CPASMA')
    ip = np.abs(ip_node.data())

    # internal inductance
    li = OMFITmdsValue(server='CMOD', shot=shot, treename='ANALYSIS',TDI='\\analysis::EFIT_AEQDSK:ali').data()

    R_cm = 67.0 # value chosen/fixed in scopes
    L = li*2.*np.pi*R_cm*1e-9  # total inductance (nH)
    
    vi = L * np.gradient(smooth(ip,2),time)   # induced voltage
    
    P_oh = ip * (vsurf - vi)/1e6 # P=IV   #MW
    return time, P_oh

    
def get_CMOD_var(var,shot, tmin=None, tmax=None, plot=False):
    ''' Get tree variable for a CMOD shot. If a time window is given, the value averaged over that window is returned,
    or else the time series is given.  See list below for acceptable input variables.
    '''

    if var=='Bt':
        node = OMFITmdsValue(server='CMOD', shot=shot, treename='magnetics',TDI='\\magnetics::Bt')
    elif var=='Bp':
        # use Bpolav, average poloidal B field --> see definition in Silvagni NF 2020
        node = OMFITmdsValue(server='CMOD',shot=shot,treename='analysis', TDI='\EFIT_AEQDSK:bpolav')
    elif var=='Ip':
        node = OMFITmdsValue(server='CMOD', shot=shot, treename='magnetics',TDI='\\magnetics::Ip')
    elif var=='nebar':
        node = OMFITmdsValue(server='CMOD', shot=shot, treename='electrons',TDI='\\electrons::top.tci.results:nl_04')
    elif var=='P_RF':
        node = OMFITmdsValue(server='CMOD', shot=shot, treename='RF',TDI='\\RF::RF_power_net')
    elif var=='P_ohmic':
        node = OMFITmdsValue(server='CMOD', shot=shot, treename='RF',TDI='\\RF::RF_power_net')
    elif var=='P_rad':
        try: 
            node = OMFITmdsValue(server='CMOD', shot=shot, treename='spectroscopy',TDI='\\spectroscopy::top.bolometer:results:foil:main_power') # W
            getdata = node.data()
            radvar = 'main'
        except:
            node = OMFITmdsValue(server='CMOD', shot=shot, treename='spectroscopy',TDI='\\spectroscopy::top.bolometer:twopi_diode') # kW
            getdata = node.data()
            radvar = 'diode'
    elif var=='p_D2':
        node = OMFITmdsValue(server='CMOD', shot=shot, treename='EDGE',TDI='\\EDGE::TOP.GAS.RATIOMATIC.F_SIDE')  # mTorr
    elif var=='p_E_BOT_MKS':
        node = OMFITmdsValue(server='CMOD', shot=shot, treename='EDGE',TDI='\\EDGE::E_BOT_MKS')  # mTorr   #lower divertor
    elif var=='p_B_BOT_MKS':
        node = OMFITmdsValue(server='CMOD', shot=shot, treename='EDGE',TDI='\\EDGE::B_BOT_MKS')     # mTorr  # lower divertor
    elif var=='p_F_CRYO_MKS':
        node = OMFITmdsValue(server='CMOD', shot=shot, treename='EDGE',TDI='\\EDGE::F_CRYO_MKS')     # mTorr, only post 2006
    elif var=='p_G_SIDE_RAT':
        node = OMFITmdsValue(server='CMOD', shot=shot, treename='EDGE',TDI='\\EDGE::G_SIDE_RAT')   
    elif var=='q95':
        node = OMFITmdsValue(server='CMOD',shot=shot, treename='analysis', TDI='\EFIT_AEQDSK:qpsib')
    elif var=='Wmhd':
        node = OMFITmdsValue(server='CMOD',shot=shot, treename='analysis', TDI='\EFIT_AEQDSK:wplasm')
    elif var=='areao':
        node = OMFITmdsValue(server='CMOD',shot=shot, treename='analysis', TDI='\EFIT_AEQDSK:areao')
    elif var=='betat':
        node = OMFITmdsValue(server='CMOD',shot=shot, treename='analysis', TDI='\EFIT_AEQDSK:betat')
    elif var=='betap':
        node = OMFITmdsValue(server='CMOD',shot=shot, treename='analysis', TDI='\EFIT_AEQDSK:betap')
    elif var=='P_oh':
        t,data = get_P_ohmic(shot)   # accurate routine to estimate Ohmic power
    elif var=='h_alpha':
        node = OMFITmdsValue(server='CMOD',shot=shot, treename='spectroscopy', TDI='\ha_2_bright')
    elif var=='cryo_on':
        node = OMFITmdsValue(server='CMOD',shot=shot, treename='EDGE', TDI='\EDGE::TOP.CRYOPUMP:MESSAGE')
    elif var=='ssep':
        node = OMFITmdsValue(server='CMOD', shot=shot, treename='analysis',TDI='\EFIT_AEQDSK:ssep')
    else:
        raise ValueError('Variable '+var+' was not recognized!')

    if var not in ['P_oh']: 
        data = node.data()
        t = node.dim_of(0)

        if var=='p_E_BOT_MKS' or var=='p_B_BOT_MKS' or var=='p_F_CRYO_MKS':  # anomalies in data storage
            data = data[0,:]
    
    if var=='P_rad':
        if radvar == 'diode': # if radvar == 'main', no need to scale
            # From B.Granetz's matlab scripts: factor of 4.5 from cross-calibration with 2pi_foil during flattop
            # NB: Bob's scripts mention that this is likely not accurate when p_rad (uncalibrated) <= 0.5 MW
            #data *= 4.5
            data *= 3 # suggestion by JWH
            # data from the twopi_diode is output in kW. Change to MW for consistency
            data /= 1e3
        if radvar == 'main': # if radvar == 'main', just need to convert to MW
            data /= 1e6

    if plot:
        plt.figure()
        plt.plot(t,data)
        plt.xlabel('time [s]')
        plt.ylabel(var)

    if tmin is not None and tmax is not None:
        tidx0 = np.argmin(np.abs(t - tmin))
        tidx1 = np.argmin(np.abs(t - tmax))
        return np.mean(data[tidx0:tidx1])
    else:
        return t,data


def load_fmp_neTe(shot, get_max=False, plot=False):
    '''Load slow ne and Te from Flush Mounted Probes (FMP) on the divertor from the nodes
    \EDGE::top.probes.fmp.osd_0{ii}.p0.*e_Slow

    If get_max=True, returns the maximum of all the loaded probe signals over time. Otherwise, return
    individual signals.
    '''
    ne_fmp = []
    Te_fmp = []
    t_fmp = []

    ii=1
    while True:
        node_ne = OMFITmdsValue(server='CMOD',shot=shot,treename='EDGE', TDI=f'\EDGE::top.probes.fmp.osd_0{ii}.p0.ne_Slow')
        node_Te = OMFITmdsValue(server='CMOD',shot=shot,treename='EDGE', TDI=f'\EDGE::top.probes.fmp.osd_0{ii}.p0.Te_Slow')

        if node_ne.data() is None:
            break
        
        ne_fmp.append(node_ne.data())
        t_fmp.append(node_ne.dim_of(0))
        Te_fmp.append(node_Te.data())
        ii+=1

    
    if get_max:
        ne_fmp_interp = np.zeros((len(ne_fmp),200)) # 200 time points is enough, usually ~100 in signals
        Te_fmp_interp = np.zeros((len(Te_fmp),200)) # 200 time points is enough, usually ~100 in signals
        
        # each probe has a different time base. Interpolate and then sum
        tmin = np.min([np.min(tlist) for tlist in t_fmp])
        tmax= np.max([np.max(tlist) for tlist in t_fmp])
        time = np.linspace(tmin, tmax, ne_fmp_interp.shape[1])
        
        for ii in np.arange(len(ne_fmp)):
            ne_fmp_interp[ii,:] = interp1d(t_fmp[ii], ne_fmp[ii], bounds_error=False)(time)
            Te_fmp_interp[ii,:] = interp1d(t_fmp[ii], Te_fmp[ii], bounds_error=False)(time)
            
        ne_fmp_max = np.nanmax(ne_fmp_interp, axis=0)
        Te_fmp_max = np.nanmax(Te_fmp_interp, axis=0)

        if plot:
            fig,ax = plt.subplots()
            ax.plot(time, ne_fmp_max)
            ax.set_xlabel('time [s]')
            ax.set_ylabel(r'$n_e$ FMP max [m$^{-3}$]')
            
            fig,ax = plt.subplots()
            ax.plot(time, Te_fmp_max)
            ax.set_xlabel('time [s]')
            ax.set_ylabel(r'$T_e$ FMP max [eV]')

        return time, ne_fmp_max, Te_fmp_max
    
    if not get_max and plot:
        fig1,ax1 = plt.subplots()
        fig2,ax2 = plt.subplots()
        
        for elem in np.arange(len(ne_fmp)):
            ax1.plot(t_fmp[elem], ne_fmp[elem], label=f'{elem}')
            ax2.plot(t_fmp[elem], Te_fmp[elem], label=f'{elem}')
            
        ax1.set_xlabel('time [s]')
        ax1.set_ylabel(r'$n_e$ [m$^{-3}$]')
        plt.tight_layout(); plt.legend()
        ax2.set_xlabel('time [s]')
        ax2.set_ylabel(r'$T_e$ [m$^{-3}$]')
        plt.tight_layout(); plt.legend()

    return t_fmp, ne_fmp, Te_fmp


def Lya_to_ion_rate(emiss_prof, ne, Te, ni=None, rhop=None, rates_source='adas'):
    '''Estimate ionization rate measured from ground state density and emissivity profiles.'''


    assert len(emiss_prof)==len(ne) and len(ne)==len(Te)
    if ni is None:
        ni=copy.deepcopy(ne)
    else:
        assert len(ne)==len(ni)

    nn, _ = aurora.Lya_to_neut_dens(emiss_prof, ne, Te, plot=False, rhop=rhop, rates_source=rates_source)

    atom_data = aurora.get_atom_data('H')
    ion_func = aurora.interp_atom_prof(atom_data['scd'], np.log10(ne), np.log10(Te), x_multiply=True)
    ion_rate = ion_func[:,0]*nn

    return ion_rate




