import sys
import numpy
import math
from matplotlib import pyplot

#galpy
from galpy.orbit import Orbit
from galpy.potential import LogarithmicHaloPotential, MWPotential2014
from galpy.util import bovy_plot, bovy_coords, bovy_conversion


# Solar position (kpc) and Local Standard of Rest (km/s)
_R0 =   8.0 # Sun is located 8 kpc away from the Galactic Center
_V0 = 220.0 # LSR velocity is 220 km/s
_z0 =   0.0 # Sun is on the Galactic mid-plane (_z0 is in kpc)

# Solar peculiar motion (km/s)
_Usun = -11.1  # _Usun < 0 means Sun is moving towards Galactic Center
_Vsun =  24.   # _Vsun > 0 means Sun is faster than Local Standard of Rest
_Wsun =   7.25 # _Wsun > 0 means Sun is moving towards North Galactic Pole

# J2000
_my_epoch = 2000.0

def example_integrate_orbit():
    # Use 'MWPotential2014' as an example
    my_potential = MWPotential2014
    
    # integration time in units of (_R0/_V0)
    time_start   =   0.
    time_end     = 100.
    my_time_step = numpy.linspace(time_start,time_end,1001)

    # read 6D astrometric data
    # (0) [mas] parallax
    # (1) [deg] ell
    # (2) [deg] b
    # (3) [km/s] heliocentric line-of-sight velocity
    # (4) [mas/yr] proper motion along ell direction (mu_ellstar = mu_ell * cos(b))
    # (5) [mas/yr] proper motion along b direction (mu_b)

    my_filename = 'parallax_ell_b_heliocentricLineOfSightVelocity_properMotionEllStar_properMotionB.txt'
    my_file     = open(my_filename, 'r')
    my_data     = numpy.loadtxt(my_file, comments='#')
    
    my_parallax_mas    = my_data[:,0]
    my_ell_deg         = my_data[:,1]
    my_b_deg           = my_data[:,2]
    my_hlosv_kms       = my_data[:,3]
    my_muellstar_masyr = my_data[:,4]
    my_mub_masyr       = my_data[:,5]
    
    
    # count sample size
    my_sample_size = len(my_ell_deg)
    
    
    for i in range (my_sample_size):
        print('star ID=%d' % (i))
        
        # convert parallax to distance
        distance_kpc = 1./my_parallax_mas[i]
    
        # convert (ell, b) to (RA, DEC)
        RA_deg, DEC_deg = bovy_coords.lb_to_radec(my_ell_deg[i], my_b_deg[i], degree=True, epoch=_my_epoch)
    
        # heliocentric line-of-sight velocity
        hlosv_kms = my_hlosv_kms[i]

        # convert (mu_ellstar, mu_b) to (mu_RAstar, mu_DEC)
        muRAstar_masyr, muDEC_masyr = bovy_coords.pmllpmbb_to_pmrapmdec(my_muellstar_masyr[i],my_mub_masyr[i], my_ell_deg[i], my_b_deg[i], degree=True, epoch=_my_epoch)
        
        # create orbit instance
        obs_6D = Orbit(vxvv=[RA_deg,DEC_deg,distance_kpc,muRAstar_masyr,muDEC_masyr,hlosv_kms],radec=True,ro=_R0,vo=_V0,zo=0.,solarmotion=[_Usun,_Vsun,_Wsun])

        # integrate orbit in my_potential
        obs_6D.integrate(my_time_step,my_potential)

        # file on which we write 6D data at each time step
        outfile_i = open("t_x_y_z_vx_vy_vz_R__orbitID%03d.txt" % (i),'w')
        
        # For illustrative purpose, I use an explicit expression to access 6D data at each time.
        for j in range (len(my_time_step)):
            # Here I assume that
            # (a) Galactic Center is located at (x,y)=(0,0) kpc
            # (b) Sun is located at (x,y)=(-8,0) kpc
            # (c) Nearby disc stars with circular orbits move towards (vx,vy)=(0,220) km/s
            # This is why I add a minus sign (-) to x and vx.
            x = -obs_6D.x (my_time_step[j])
            y =  obs_6D.y (my_time_step[j])
            z =  obs_6D.z (my_time_step[j])
            vx= -obs_6D.vx(my_time_step[j])
            vy=  obs_6D.vy(my_time_step[j])
            vz=  obs_6D.vz(my_time_step[j])
            R =  obs_6D.R (my_time_step[j])
            printline = '%lf %lf %lf %lf %lf %lf %lf %lf\n' % (my_time_step[j],x,y,z,vx,vy,vz,R)
            outfile_i.write(printline)

        # close file
        outfile_i.close()

    return None

if __name__ == '__main__':
    example_integrate_orbit()
