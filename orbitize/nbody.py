import os
import numpy as np
import astropy
import astropy.units as u
import orbitize.basis as basis
import rebound
import matplotlib.pylab as plt

# TODO: add optional keyword argument that is an array of masses of planets (m_pl)
# TODO: if m_pl is None, then assume masses of planets are zero
# TODO: otherwise, m_star = mtot - sum(m_pl)

def calc_orbit(epochs, sma, ecc, inc, aop, pan, tau, plx, mtot, tau_ref_epoch):
    """
    Solves for position for a set of input orbital elements using rebound.

    Args:
        epochs (np.array): MJD times for which we want the positions of the planet
        sma (np.array): semi-major axis of orbit [au]
        ecc (np.array): eccentricity of the orbit [0,1]
        inc (np.array): inclination [radians]
        aop (np.array): argument of periastron [radians]
        pan (np.array): longitude of the ascending node [radians]
        tau (np.array): epoch of periastron passage in fraction of orbital period past MJD=0 [0,1]
        plx (np.array): parallax [mas]
        mtot (np.array): total mass of the two-body orbit (M_* + M_planet) [Solar masses]
        tau_ref_epoch (float, optional): reference date that tau is defined with respect to (i.e., tau=0)

    Returns:
        3-tuple:

            raoff (np.array): array-like (n_dates x n_orbs) of RA offsets between the bodies
            (origin is at the other body) [mas]

            deoff (np.array): array-like (n_dates x n_orbs) of Dec offsets between the bodies [mas]

            vz (np.array): array-like (n_dates x n_orbs) of radial velocity of one of the bodies
                (see `mass_for_Kamp` description)  [km/s]
        return is in format [raoff[planet1, planet2,...,planetn], deoff[planet1, planet2,...,planetn], vz[planet1, planet2,...,planetn]
    """

    # initialze a 2-body system with the input orbital parameters, starting at the first epoch
    # run the simulation forward in time until the last epoch
    # return the position & velocity of the planet at each of the epochs

    
    sim = rebound.Simulation()      #creating the simulation in Rebound
    sim.units = ('AU','yr','Msun')  #converting units for uniformity
    sim.add(m = mtot)               #give all mass to star, planet mass = 0, for now
    ps = sim.particles              #for easier calls

    tx = len(epochs)                #keeping track of how many time steps
    te = epochs-epochs[0]

    indv = len(sma)                 #number of planets orbiting the star
    num_planets = np.arange(0,indv) #creates an array of indeces for each planet that exists
    ra_reb = np.zeros((tx, indv))   #numpy.zeros(number of [arrays], size of each array)
    dec_reb = np.zeros((tx, indv))
    vz = np.zeros((tx, indv))

    #for each planet, create a body in the Rebound sim
    for i in num_planets:
        mnm = basis.tau_to_manom(epochs[0], sma[i], mtot, tau[i], tau_ref_epoch) #calculating mean anomaly
        sim.add(m = 0, a = sma[i], e = ecc[i], inc = inc[i], Omega = pan[i] + np.pi/2, omega =aop[i], M =mnm)
        sim.move_to_com()
        sim.integrator = "ias15" 
        sim.dt = ps[1].P/1000.

    #integrate at each epoch
    for j,t in enumerate(te):
        sim.integrate(t/365.25)

        #for each planet in each epoch denoted by j,t find the RA, Dec, and RV
        for i in num_planets:
            ra_reb[j,i] = -(ps[int(i+1)].x - ps[0].x) # ra is negative x
            dec_reb[j,i] = ps[int(i+1)].y - ps[0].y
            vz[j,i] = ps[int(i+1)].vz

    #adjusting for parallax
    pxr = plx*ra_reb 
    pxd = plx*dec_reb

    if len(sma)==1:                 #for formatting purposes
        planetra = pxr.reshape(tx,)
        planetdec = pxd.reshape(tx,)
        planetvz = vz.reshape(tx,)
        return planetra, planetdec, planetvz

    else: 
        return pxr, pxd, vz

