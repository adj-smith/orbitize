

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

    
    import os
    import numpy as np
    import astropy.table as table
    import astropy
    import astropy.units as u
    import orbitize
    import orbitize.read_input as read_input
    import orbitize.kepler as kepler
    import orbitize.system as system
    import orbitize.basis as basis
    import orbitize.priors as priors
    import orbitize.driver as driver
    import rebound

"""
For Sofia:

orbitize.kepler.calc_orbit(epochs, sma, ecc, inc, aop, pan, tau, plx, mtot, 
    mass_for_Kamp=None, tau_ref_epoch=0, tolerance=1e-09, max_iter=100)

rebound.OrbitPlot(sim, figsize=None, fancy=False, slices=0, xlim=None, ylim=None, 
    unitlabel=None, color=False, periastron=False, orbit_type='trail', lw=1.0, 
    plotparticles=[], primary=None, Narc=128)

"""
#add additional arg for mstar/mplanet?
sim = rebound.Simulation()
sim.units = ('AU','days','Msun')
sim.add(m = mtot-somem )
sim.add(m = , a = , e = , inc = , Omega = , omega = )






#for temp. testing purposes only
def test_1planet():
    """
    Sanity check that things agree for 1 planet case
    """
    # generate a planet orbit
    sma = 1
    ecc = 0.1
    inc = np.radians(45)
    aop = np.radians(45)
    pan = np.radians(45)
    tau = 0.5
    plx = 1
    mtot = 1
    tau_ref_epoch = 0
    mjup = u.Mjup.to(u.Msun)
    mass_b = 12 * mjup

    epochs = np.linspace(0, 300, 100) + tau_ref_epoch # nearly the full period, MJD

    ra_model, dec_model, vz_model = kepler.calc_orbit(epochs, sma, ecc, inc, aop, pan, tau, plx, mtot, tau_ref_epoch=tau_ref_epoch)

    # generate some fake measurements just to feed into system.py to test bookkeeping
    t = table.Table([epochs, np.ones(epochs.shape, dtype=int), ra_model, np.zeros(ra_model.shape), dec_model, np.zeros(dec_model.shape)], 
                     names=["epoch", "object" ,"raoff", "raoff_err","decoff","decoff_err"])
    filename = os.path.join(orbitize.DATADIR, "rebound_1planet.csv")
    t.write(filename)

    # create the orbitize system and generate model predictions using the ground truth
    astrom_dat = read_input.read_file(filename)

    sys = system.System(1, astrom_dat, mtot, plx, tau_ref_epoch=tau_ref_epoch)

    params = np.array([sma, ecc, inc, aop, pan, tau, plx, mtot])
    radec_orbitize, _ = sys.compute_model(params)
    ra_orb = radec_orbitize[:, 0]
    dec_orb = radec_orbitize[:, 1]


    # now project the orbit with rebound
    manom = basis.tau_to_manom(epochs[0], sma, mtot, tau, tau_ref_epoch)
    
    sim = rebound.Simulation()
    sim.units = ('yr', 'AU', 'Msun')


    # add star
    sim.add(m=mtot - mass_b)

    # add one planet
    sim.add(m=mass_b, a=sma, e=ecc, M=manom, omega=aop, Omega=pan+np.pi/2, inc=inc)
    ps = sim.particles

    sim.move_to_com()

    # Use Wisdom Holman integrator (fast), with the timestep being < 5% of inner planet's orbital period
    sim.integrator = "ias15"
    sim.dt = ps[1].P/1000.

    # integrate and measure star/planet separation 
    ra_reb = []
    dec_reb = []

    for t in epochs:
        sim.integrate(t/365.25)
        
        ra_reb.append(-(ps[1].x - ps[0].x)) # ra is negative x
        dec_reb.append(ps[1].y - ps[0].y)
        
    ra_reb = np.array(ra_reb)
    dec_reb = np.array(dec_reb)

    diff_ra = ra_reb - ra_orb/plx
    diff_dec = dec_reb - dec_orb/plx

    assert np.all(np.abs(diff_ra) < 1e-9)
    assert np.all(np.abs(diff_dec) < 1e-9)