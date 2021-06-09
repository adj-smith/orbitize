import os
import numpy as np
import astropy
import astropy.units as u
import orbitize.basis as basis
import rebound
import matplotlib.pylab as plt

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

    
    sim = rebound.Simulation() #creating the simulation in Rebound
    sim.units = ('AU','yr','Msun') 
    sim.add(m = mtot) #give all mass to star, planet mass = 0, for now
    ps = sim.particles #for easier calls

    tx = len(epochs)
    te = epochs-epochs[0]   

    indv = len(sma)
    ra_reb = np.zeros((tx, indv))
    dec_reb = np.zeros((tx, indv))
    vz = np.zeros((tx, indv))

    for i in np.arange(0,indv):
        mnm = basis.tau_to_manom(epochs[0], sma[i], mtot, tau[i], tau_ref_epoch) #calculating mean anomaly
        sim.add(m = 0, a = sma[i], e = ecc[i], inc = inc[i], Omega = pan[i] + np.pi/2, omega =aop[i], M =mnm)
        
        sim.move_to_com()
        sim.integrator = "ias15" 
        sim.dt = ps[1].P/1000.

        for j,t in enumerate(te):
            #print(j,t)
            print(type(i+1))
            sim.integrate(t/365.25)
            ra_reb[j,i] = -(ps[int(i+1)].x - ps[0].x) # ra is negative x
            dec_reb[j,i] = ps[int(i+1)].y - ps[0].y
            vz[j,i] = ps[int(i+1)].vz

        
 

    pxr = plx*ra_reb #adjusting for parallax
    pxd = plx*dec_reb #adjusting for parallax

    return pxr, pxd, vz

#Test Data
"""
sma = np.array([1])
ecc = np.array([0.1])
inc = np.array([np.radians(45)])
aop = np.array([np.radians(45)])
pan = np.array([np.radians(45)])
tau = np.array([0.5])
plx = np.array([1])
mtot = np.array([1])
tau_ref_epoch = np.array([0])
epochs = np.linspace(0, 1000, 1000) + tau_ref_epoch # nearly the full period, MJD
"""

import orbitize.kepler

#From System HR-8799
epochs = np.linspace(0, 1000, 1000)
sma = np.array([71, 43, 26, 16])
ecc = np.array([.02, .02, .13, .12])
inc = np.array([0.47123, 0.47123, 0.47123, 0.47123])
aop = np.array([1.5184, 1.16937, 0.29670, 1.91986])
pan = np.array([1.18682, 1.18682, 1.18682, 1.18682])
tau = np.array([0.54, 0.50, 0.79, 0.71])
plx = np.array([7])
mtot = np.array([1.49])
tau_ref_epoch = 0

"""
x1,x2,x3 = (calc_orbit(epochs,sma,ecc,inc,aop,pan,tau,plx,mtot, tau_ref_epoch))
y1,y2,y3 = (orbitize.kepler.calc_orbit(epochs,sma,ecc,inc,aop,pan,tau,plx,mtot, tau_ref_epoch))
print('rebound:', x1[:7],',', x2[57],',',x3[57])
print('Kepler:', y1[:7],',', y2[57],',', y3[57])
"""

#Analysis Functions

def time_analysis(nt):
    import time
    import orbitize.kepler
    import matplotlib.pyplot as plt

    #where nt is the longest time you want to run the solver for, in epochs
    te = np.linspace(0,nt,1000) #range of number of epochs you want to test
    index_name = np.arange(0,len(te)) #simple index for following 'for' loop

    treb = np.zeros(len(te)) #time rebound takes, loadable array
    tkep = np.zeros(len(te)) #time kepler takes, loadable array
    num_epochs = np.zeros(len(te)) #number of epochs measured in each calculation, will be the x axis

    for i in index_name:
        var_epoch = np.linspace(0, te[i], 100) + tau_ref_epoch #what will be plugged into each solver
        num_e = te[i]/100 #number of epochs used
        num_epochs[i] = num_e #number of epochs for each point

        t1 = time.time()
        calc_orbit(var_epoch, sma,ecc,inc,aop,pan,tau,plx,mtot,tau_ref_epoch)
        t2 = time.time()
        orbitize.kepler.calc_orbit(var_epoch, sma,ecc,inc,aop,pan,tau,plx,mtot, tau_ref_epoch)
        t3 = time.time()
        tr = t2-t1
        tk = t3-t2

        treb[i] = tr
        tkep[i] = tk

    plt.plot(num_epochs, treb, label = 'Rebound')
    plt.plot(num_epochs, tkep, label = 'Current Keplerian')
    plt.ylabel('Time it takes to solve (seconds)')
    plt.xlabel('Number of Epochs')

    axes1 = plt.gca()
    axes2 = axes1.twiny()
    axes2.set_xticks([.33, .66, .99])
    axes1.set_xlabel("Number of Epochs")
    axes2.set_xlabel("Total time calculated (Earth years)")
 

    plt.legend()
    plt.show()


def calc_diff():
    import orbitize.kepler
    import matplotlib.pyplot as plt

    rra, rde, rvz = (calc_orbit(epochs, sma,ecc,inc,aop,pan,tau,plx,mtot,tau_ref_epoch))
    kra, kde, kvz = (orbitize.kepler.calc_orbit(epochs, sma,ecc,inc,aop,pan,tau,plx,mtot, tau_ref_epoch))

    delta_ra = abs(rra-kra)
    delta_de = abs(rde-kde)
    delta_vz = abs(rvz-kvz)

    """
    print(rvz)
    print(kvz)
    """
    Ar = np.zeros(1000) #planet A
    Ad = np.zeros(1000)
    Av = np.zeros(1000)

    Br = np.zeros(1000) #planet B
    Bd = np.zeros(1000)
    Bv = np.zeros(1000)
    
    Cr = np.zeros(1000) #planet C
    Cd = np.zeros(1000)
    Cv = np.zeros(1000)
    
    Dr = np.zeros(1000) #planet D
    Dd = np.zeros(1000)
    Dv = np.zeros(1000)


    for i in np.arange(0,999):
        i2 = int(i)
        
        Ar[i] = delta_ra[i2,0] #first planet
        Ad[i] = delta_de[i2,0]
        Av[i] = delta_vz[i2,0]
        
        Br[i] = delta_ra[i2,1] #second planet
        Bd[i] = delta_de[i2,1]
        Bv[i] = delta_vz[i2,1]

        Cr[i] = delta_ra[i2,2] #third planet
        Cd[i] = delta_de[i2,2]
        Cv[i] = delta_vz[i2,2]

        Dr[i] = delta_ra[i2,3] #fourth planet
        Dd[i] = delta_de[i2,3]
        Dv[i] = delta_vz[i2,3]
        

    plt.plot(epochs, Ar, 'brown', label = 'Planet A: RA offsets') #first planet
    plt.plot(epochs, Ad, 'red', label = 'Planet A: Dec offsets')
    #plt.plot(epochs, Av, 'pink', label = 'Planet A: RV offsets')

    plt.plot(epochs, Br, 'coral', label = 'Planet B: RA offsets') #second planet
    plt.plot(epochs, Bd, 'darkorange', label = 'Planet B: Dec offsets')
    #plt.plot(epochs, Bv, 'gold', label = 'Planet B: RV offsets')

    plt.plot(epochs, Cr, 'greenyellow', label = 'Planet C: RA offsets') #third planet
    plt.plot(epochs, Cd, 'green', label = 'Planet C: Dec offsets')
    #plt.plot(epochs, Cv, 'darkgreen', label = 'Planet C: RV offsets')

    plt.plot(epochs, Dr, 'dodgerblue', label = 'Planet D: RA offsets') #fourth planet
    plt.plot(epochs, Dd, 'blue', label = 'Planet D: Dec offsets')
    #plt.plot(epochs, Dv, 'indigo', label = 'Planet D: RV offsets')

    plt.xlabel('Epochs (Earth years)')
    plt.ylabel('Difference between Kepler and N-Body solver (milliarcseconds)')
    plt.legend()
    plt.show()
    
    
