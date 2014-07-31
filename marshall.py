#! /usr/bin/env python

from uncertainties import ufloat
import matplotlib.pyplot as plt
import numpy as np

import csv
import argparse

__FILENAME__ = "marshall.dat"

def q_marshall(lon, lat):
    """
        Looks up lat and long in the marshall file and returns an array of the
        radii and of the extinctions at each radii.

        Parameters
        ---------
                lat : float
                Latitude of object to lookup

                lon : float
                Longitude of object to lookup

        Returns
        ---------
                ext : list
                List of distance dependent extinction values and error i.e. a 
                list of ufloats

                rad : list
                List of radii of the extinction cuts with uncertainty i.e. a 
                list of ufloats
    """

    lon, lat = conCoords(lon, lat)

    with open(__FILENAME__, 'r') as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            pos = [float(row[0]), float(row[1])]
            if (lon == pos[0]) & (lat == pos[1]):
                ext, rad = split_dat(row[3:])

    f.close()

    return ext, rad

def conCoords(lon, lat):
    """
        Converts the longitude and latitude into the correct format.

        Parameters
        ----------
                lon : float
                The longitude in degrees.

                lat : float
                The latitude in degrees.

        Returns
        ----------
                lon : float
                Longitude in the correct format. -100 deg < l < 100 deg.

                lat : float
                Latitude in the correct format. |b| < 10 deg.
    """
    lon = (lon - 360) if lon > 180.0 else lon

    assert (abs(lon) <= 100) and (abs(lat) <= 10), 'Valid latitude and longitude range; |b| <= 10 deg, |l| <= 100 deg'

# if gets here then lon and lat in form and range -100 < l < 100 and |b| < 10 deg

    lon = round(lon*4)/4 if (lon % 0.25 != 0) else lon

    lon = 360 + lon if (lon < 0) else lon

    lat = round(lat*4)/4  if (lat % 0.25 != 0) else lat

    return lon, lat

def split_dat(data):
    """
        Takes the row in the form r_n e_r_n ext_n e_ext_n and returns two lists
        one with the radii and the other with the extinction values both lists
        are of ufloats

        Parameters
        --------
                data : list
                A list of floats of the radii and errors and extinctions and 
                errors.

        Returns 
        --------
                rad : list
                List of ufloats of the radii.

                ext : list
                List of ufloats of the extinctions.
    """
    rad, ext = [], []

    for i in range(0, len(data) - 1, 4):
        vec = [data[i], data[i+1], data[i+2], data[i+3]]
        
        if vec[0].strip() != '':
            radius = ufloat(vec[0], vec[1])
            extinct = ufloat(vec[2], vec[3])

            rad.append(radius)
            ext.append(extinct)

    return ext, rad

def find_asymptotic_red(ext, rad):
    """
        Finds the asymptotic line of sight reddening for a given set
        of extinctions and radial slices.

        Parameters
        ---------
                ext : list of ufloats
                A list of the extinction values in the k_s-band.

                rad : list of ufloats
                The list of radial slices where the extinction has been determined.

        Returns
        ---------
                asymptotic_red : ufloat
                The asymptotic line of sight reddening if the reddening 
                does not asymptote then returns None.
    """

    ext = [x / 0.114 for x in ext]             # v-band extinction

    slopes = []

    for i in range(len(ext)-1):
        slopes.append((ext[i+1] - ext[i])/ (rad[i+1] - rad[i]))

    for i in range(len(slopes)):
        if slopes[i] <= 0.01:
            return ext[i], rad[i]
            
    return None, None

def plot_ext(ext, rad, lat, lon):
    """
        Plots the extinction as a function of radius.

        Parameters
        ---------
                ext : list
                A list of ufloats of the extinction values.

                rad : list
                A list of ufloats of the radius values.

        Returns
        ---------
                None
    """

    fig, ax = plt.subplots(nrows=1, ncols=2)

    info = r"$l={}^\circ ,\, b={}^\circ$".format(lon,lat)

    ax[0].plot(1,1)

    ax[0].errorbar([0.] + [x.n for x in rad], [0.] + [x.n for x in ext], yerr=[0.] + [x.s for x in ext], xerr=[0.] + [x.s for x in rad], marker='o', ms=5, ecolor='red', label=info)

    ax[0].set_ylim(0, max([x.n for x in ext]) * 1.5)
    ax[0].set_xlim(0, max([x.n for x in rad]) * 1.5)
    ax[0].set_ylabel(r'$A_{k_s}\, {\rm [mag]}$', fontsize=16)
    ax[0].set_xlabel(r'${\rm Distance\, [kpc]}$', fontsize=16)
    ax[0].legend(loc='best', numpoints=1, fancybox=False, shadow=False)

    ext = [x / 0.114 for x in ext]

    ax[1].errorbar([0.] + [x.n for x in rad], [0.] + [x.n for x in ext], yerr=[0.] + [x.s for x in ext], xerr=[0.] + [x.s for x in rad], marker='o', ms=5, ecolor='red', label=info)

    ax[1].set_xlim(0, max([x.n for x in rad]) * 1.5)
    ax[1].set_ylabel(r'$A_{V}\, {\rm [mag]}$', fontsize=16)
    ax[1].set_xlabel(r'${\rm Distance\, [kpc]}$', fontsize=16)
    ax[1].legend(loc='best', numpoints=1, fancybox=False, shadow=False)

    plt.show()

def get_ext(lon, lat):
    """
        Function used to query the Marshall et al. (2006) 3D interstellar
        extinction distribution.

        Parameters
        ----------
                lon : float
                The longitude to query in degrees.

                lat : float
                The latitude to query in degrees.

        Returns
        ---------
                None
    """

    ext, rad = q_marshall(lon, lat)

    asymp_red, asymp_dist = find_asymptotic_red(ext, rad)

    print( "3D Galactic Reddening Map of Marshall et al. (2006)" )
    print( "Reddening determined in the K_s band from Marshall+" )
    print( "Then converted to visual extinction using CCM89 using;" )
    print( "\t\tA_ks = 0.114 x A_v" )
    print( "----------------------------------------------------------" )
    print( "Finding Extinction at; l = {} deg, b = {} deg".format(lon, lat) )
    if (not asymp_red) or (not asymp_dist):
        print( "No Asymptotic Reddening determine, possibly due to there being too few points." )
    else:
        z_height = np.abs(asymp_dist * np.sin((np.pi / 180) * lat) * 1000)
        print( "Asymptotic Reddening; A_V = {:.2F} mag @ a distance of {:.1F} kpc".format(asymp_red, asymp_dist) )
        print( "Z-height of Asymptotic Reddening; z = {::0F} pc".format(z_height) )
    print( "----------------------------------------------------------" )

    plot_ext(ext, rad, lat, lon)


def main():
    """
        Function used to query the Marshall et al. (2006) 3D interstellar
        extinction distribution, when the script is run directly from the
        command line. Input parameters are retrieved from command line
        arguments.

        Parameters
        ----------
                None

        Returns
        ---------
                None
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("-l", action='store', dest='lon', type=float, help='Longitude in degrees')
    parser.add_argument("-b", action='store', dest='lat', type=float, help='Latitude in degrees')

    args = parser.parse_args()

    if (not args.lon) or (not args.lat):
        print( 'Latitude and Longitude not defined' )
        quit()
    else:
        lon = args.lon
        lat = args.lat

    get_ext(lon, lat)

if __name__ == '__main__':
    main()
