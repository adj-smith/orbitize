import os
import numpy as np
import astropy
import astropy.units as u
import orbitize.basis as basis
import orbitize.kepler
import rebound
import matplotlib.pyplot as plt

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
    sim.G = 39.476926408897626                   #
    sim.add(m = mtot)               #give all mass to star, planet mass = 0, for now
    ps = sim.particles              #for easier calls

    tx = len(epochs)                #keeping track of how many time steps
    te = epochs-epochs[0]

    indv = len(sma)                 #number of planets orbiting the star
    num_planets = np.arange(0,indv) #creates an array of indeces for each planet that exists
    ra_reb = np.zeros((tx, indv))   #numpy.zeros(number of [arrays], size of each array)
    dec_reb = np.zeros((tx, indv))
    vz = np.zeros((tx, indv))

    P = np.zeros(len(sma))

    #for each planet, create a body in the Rebound sim
    for i in num_planets:
        mnm = basis.tau_to_manom(epochs[0], sma[i], mtot, tau[i], tau_ref_epoch) #calculating mean anomaly
        sim.add(m = 0, a = sma[i], e = ecc[i], inc = inc[i], Omega = pan[i] + np.pi/2, omega =aop[i], M =mnm)
        sim.move_to_com()
        sim.integrator = "ias15" 
        sim.dt = ps[1].P/1000. #good rule of thumb: timestep should be at most 10% of the shortest orbital period

        P[i] = ps[int(i+1)].P

    #integrate at each epoch
    for j,t in enumerate(te):
        sim.integrate(t/365.25)

        #for each planet in each epoch denoted by j,t find the RA, Dec, and RV
        for i in num_planets:
            ra_reb[j,i] = -(ps[int(i+1)].x - ps[0].x) # ra is negative x
            dec_reb[j,i] = ps[int(i+1)].y - ps[0].y
            vz[j,i] = ps[int(i+1)].vz

    #adjusting for parallax
    raoff = plx*ra_reb 
    deoff = plx*dec_reb

    #for formatting purposes
    if len(sma)==1:
        raoff = raoff.reshape(tx,)
        deoff = deoff.reshape(tx,)
        vz = vz.reshape(tx,)
        return raoff, deoff, vz, P

    else: 
        return raoff, deoff, vz, P

#Test Data
#From System HR-8799
#NOTE planets get closer to the star in alphabetical order, i.e. B is farthest, E is closest
sma = np.array([26]) 
ecc = np.array([.13])
inc = np.array([0.47123])
aop = np.array([0.29670])
pan = np.array([1.18682])
tau = np.array([0.79])
plx = np.array([7])
mtot = np.array([1.49])
tau_ref_epoch = 0
years = 365.25*300
epochs = np.linspace(0,years,1000)

import orbitize.kepler

def calc_diff():
    import orbitize.kepler
    import matplotlib.pyplot as plt

    rra, rde, rvz, P = calc_orbit(epochs, sma,ecc,inc,aop,pan,tau,plx,mtot,tau_ref_epoch)
    kra, kde, kvz = orbitize.kepler.calc_orbit(epochs, sma,ecc,inc,aop,pan,tau,plx,mtot, tau_ref_epoch)
        
    delta_ra = -(rra-kra)
    delta_de = -(rde-kde)
    #delta_vz = abs(rvz-kvz)
    yepochs = epochs/365.25

    if len(sma)==1:
        plt.plot(yepochs, delta_ra, label = 'Planet B: RA offsets')
        plt.plot(yepochs, delta_de, label = 'Planet B: Dec offsets')
        #plt.plot(yepochs, delta_vz, 'pink', label = 'Planet B: RV offsets')

    elif len(sma)==4:

        plt.plot(yepochs, delta_ra[:,0], 'brown', label = 'Planet B: RA offsets') #first planet
        plt.plot(yepochs, delta_de[:,0], 'red', label = 'Planet B: Dec offsets')
        #plt.plot(yepochs, delta_vz[:,0], 'pink', label = 'Planet B: RV offsets')

        plt.plot(yepochs, delta_ra[:,1], 'coral', label = 'Planet C: RA offsets') #second planet
        plt.plot(yepochs, delta_de[:,1], 'orange', label = 'Planet C: Dec offsets')
        #plt.plot(yepochs, delta_vz[:,1], 'gold', label = 'Planet C: RV offsets')

        plt.plot(yepochs, delta_ra[:,2], 'greenyellow', label = 'Planet D: RA offsets') #third planet
        plt.plot(yepochs, delta_de[:,2], 'green', label = 'Planet D: Dec offsets')
        #plt.plot(yepochs, delta_vz[:,2], 'darkgreen', label = 'Planet D: RV offsets')            

        plt.plot(yepochs, delta_ra[:,3], 'dodgerblue', label = 'Planet E: RA offsets') #fourth planet
        plt.plot(yepochs, delta_de[:,3], 'blue', label = 'Planet E: Dec offsets')
        #plt.plot(yepochs, delta_vz[:,3], 'indigo', label = 'Planet E: RV offsets')

    else:
        print('I dont feel like it')



    plt.xlabel('Epochs (Earth years)')
    plt.ylabel('Difference between Kepler and N-Body solver (milliarcseconds)')
    plt.legend()
    plt.show()   

    
def P_test(x):
    """
    Args: 
        x = 0: Plots RA from both Kepler and Rebound over time 
        x = 1: plots delta RA over time

    """
    rra, rde, rvz, P = calc_orbit(epochs, sma,ecc,inc,aop,pan,tau,plx,mtot,tau_ref_epoch)
    kra, kde, kvz = orbitize.kepler.calc_orbit(epochs, sma,ecc,inc,aop,pan,tau,plx,mtot, tau_ref_epoch)
    
    yepochs = epochs/365.25
    croissant = ((yepochs[-1])//P)+1
    butter = np.zeros(int(croissant))
    for i in np.arange(0,croissant):
        butter[int(i)] = i*P
    
    
    if x == 0:
        plt.plot(yepochs, rra, label = 'Planet D: Rebound', alpha = 0.5)
        plt.plot(yepochs, kra, label = 'Planet D: Kepler', alpha = 0.5)
        plt.title('Right Ascencion Offsets')
        plt.ylabel('milliarcseconds')
        plt.hlines(rra[0],yepochs[0],yepochs[-1],'gainsboro', label='RA at Epochs(0)')
        plt.vlines(butter,rra.min(), rra.max(),'gray', label='Calculated periods')


    elif x == 1:
        rra2, rde2, rvz2, P2 = calc_orbit(butter, sma,ecc,inc,aop,pan,tau,plx,mtot,tau_ref_epoch)
        kra2, kde2, kvz2 = orbitize.kepler.calc_orbit(butter, sma,ecc,inc,aop,pan,tau,plx,mtot, tau_ref_epoch)
        bread = len(rra2)
        REBdelta_ra = np.zeros(bread)
        REBdelta_de = np.zeros(bread)
        KEPdelta_ra = np.zeros(bread)
        KEPdelta_de = np.zeros(bread)

        for i in np.arange(0,bread):
            REBdelta_ra[i] = rra[i]-rra[0]
            REBdelta_de[i] = rde[i]-rde[0]
            KEPdelta_ra[i] = kra[i]-kra[0]
            KEPdelta_de[i] = kde[i]-kde[0]

        plt.plot(butter, REBdelta_ra, 'brown', label = 'Rebound delta RA offsets', alpha = 0.5)
        plt.plot(butter, REBdelta_de, 'red', label = 'Rebound delta Dec offsets', alpha = 0.5)

        plt.plot(butter, KEPdelta_ra, 'dodgerblue', label = 'Kepler delta RA offsets', alpha = 0.5)
        plt.plot(butter, KEPdelta_de, 'blue', label = 'Kepler delta Dec offsets', alpha = 0.5)

        plt.title('HR-8799D')
        plt.ylabel('Difference from initial position(milliarcseconds)')
        plt.vlines(butter,-8,2,'gray', label='Calculated periods')

    else: 
        print('Please follow instructions')

    plt.xlabel('Epochs (Earth years)')
    plt.legend()
    plt.show()

calc_diff()