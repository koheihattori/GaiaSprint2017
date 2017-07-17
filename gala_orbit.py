# Standard library
import sys
import math

# Third-party
import matplotlib.pyplot as plt
import numpy as np

# astropy and Gala
import astropy.coordinates as coord
import astropy.units as u
from gala.potential import MilkyWayPotential
from gala.dynamics import PhaseSpacePosition
from gala.units import galactic

# vx,vy,vz (or UVW), U positive *towards* the galactic center
v_sun = [11.1, 244, 7.24] * u.km/u.s

galcen_frame = coord.Galactocentric(
    galcen_distance=8*u.kpc, # Sun-Galactic center distance
    galcen_v_sun=coord.CartesianDifferential(v_sun),
    z_sun=0*u.pc) # height of the sun above the midplane

def example_integrate_orbit():
    # Use 'MilkyWayPotential' as an example
    potential = MilkyWayPotential()

    # read 6D astrometric data
    # (0) [mas] parallax
    # (1) [deg] ell
    # (2) [deg] b
    # (3) [km/s] heliocentric line-of-sight velocity
    # (4) [mas/yr] proper motion along ell direction (mu_ellstar = mu_ell * cos(b))
    # (5) [mas/yr] proper motion along b direction (mu_b)

    filename = 'parallax_ell_b_heliocentricLineOfSightVelocity_properMotionEllStar_properMotionB.txt'
    data = np.genfromtxt(filename)

    gal = coord.Galactic(distance=1000. / data[:,0] * u.pc, # parallax to distance
                         l=data[:,1] * u.deg,
                         b=data[:,2] * u.deg,
                         radial_velocity=data[:,3] * u.km/u.s,
                         pm_l_cosb=data[:,4] * u.mas/u.yr,
                         pm_b=data[:,5] * u.mas/u.yr)

    galcen = gal.transform_to(galcen_frame)
    w0 = PhaseSpacePosition(galcen.data)

    orbit = potential.integrate_orbit(w0,
                                      dt=0.5*u.Myr, # timestep
                                      n_steps=10000) # run for 2.5 Gyr
    # can also do:
    # orbit = potential.integrate_orbit(w0,
    #                                   t1=0.*u.Myr,
    #                                   t2=2.5*u.Gyr, # run for 2.5 Gyr
    #                                   n_steps=10000)

    # the 'orbit' object actually contains 2 orbits, since we specified 2 sets
    # of coordinates (from the provided text file). We can access them
    # individually by doing, e.g.:
    orbit1 = orbit[0]

    # You can now do a lot with the Orbit object that is returned. By default,
    # it is represented in Cartesian coordinates. But we can change to other
    # representations, e.g., cylindrical:
    cyl_orbit = orbit.represent_as('cylindrical')
    print(cyl_orbit.rho, cyl_orbit.v_rho)

    # Or, to plot the orbits in the Meridional plane:
    cyl_orbit.plot(['rho', 'z'])
    plt.show()

if __name__ == '__main__':
    example_integrate_orbit()
