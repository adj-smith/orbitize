import os
import numpy as np
import astropy
import astropy.units as u
import orbitize.basis as basis
import rebound
import matplotlib.pylab as plt

def calc_orbit(epochs, sma, ecc, inc, aop, pan, tau, plx, mtot, tau_ref_epoch=0):
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

    Returns:
        3-tuple:

            raoff (np.array): array-like (n_dates x n_orbs) of RA offsets between the bodies
            (origin is at the other body) [mas]

            deoff (np.array): array-like (n_dates x n_orbs) of Dec offsets between the bodies [mas]

            vz (np.array): array-like (n_dates x n_orbs) of radial velocity of one of the bodies
                (see `mass_for_Kamp` description)  [km/s]

    """

    # initialze a 2-body system with the input orbital parameters, starting at the first epoch
    # run the simulation forward in time until the last epoch
    # return the position & velocity of the planet at each of the epochs
    

    mnm = basis.tau_to_manom(epochs[0], sma, mtot, tau, tau_ref_epoch)
    sim = rebound.Simulation()
    sim.units = ('AU','yr','Msun')
    sim.add(m = mtot) #give all mass to star, planet mass = 0
    sim.add(m = 0, a = sma, e = ecc, inc = inc, Omega = pan + np.pi/2, omega =aop, M =mnm)
    sim.move_to_com()

    ps = sim.particles #for easier calls

    sim.integrator = "ias15" 
    sim.dt = ps[1].P/1000.

    tx = len(epochs)
    te = epochs-epochs[0]
    ra_reb = np.zeros(tx)
    dec_reb = np.zeros(tx)
    vz = np.zeros(tx)

    for i,t in enumerate(te):
        print(i,t)
        sim.integrate(t/365.25)
        
        ra_reb[i] = -(ps[1].x - ps[0].x) # ra is negative x
        dec_reb[i] = ps[1].y - ps[0].y
        vz[i] = ps[1].vz


    return plx*ra_reb, plx*dec_reb, vz


sma = 1
ecc = 0.1
inc = np.radians(45)
aop = np.radians(45)
pan = np.radians(45)
tau = 0.5
plx = 1
mtot = 1
tau_ref_epoch = 0
epochs = np.linspace(0, 300, 5) + tau_ref_epoch # nearly the full period, MJD

import orbitize.kepler

print('rebound test: ',calc_orbit(epochs, sma,ecc,inc,aop,pan,tau,plx,mtot,tau_ref_epoch))
print('Kepler: ',orbitize.kepler.calc_orbit(epochs, sma,ecc,inc,aop,pan,tau,plx,mtot,tau_ref_epoch))

   



    