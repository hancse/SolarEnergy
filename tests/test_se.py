
# import context  # to use local solarenergy package

import numpy as np
import pandas as pd
import pytz
import solarenergy as se


from pvlib.solarposition import get_solarposition, nrel_earthsun_distance
from pvlib.atmosphere import get_relative_airmass, alt2pres, get_absolute_airmass
from pvlib.irradiance import get_extra_radiation

import matplotlib.pyplot as plt
import matplotlib

matplotlib.use("Qt5Agg")


def test_extinction():
    am = np.linspace(1.0, 40.0, 40)
    ef = se.extinction_factor(am)
    fig, ax = plt.subplots(1)
    ax.plot(am, ef)
    ax.set_xlim(0, 40.1)
    ax.set_xlabel('Air mass / AM0')
    ax.set_ylabel('Extinction factor')
    ax.axvline(x=38.2, color='r')
    plt.title('Extinction factor vs. air mass')
    plt.show()


def test_positions():
    """VALIDATION of solar position calculations in solarenergy
       comparison against PVLIB get_solarposition()
    """
    nl_tz = pytz.timezone("Europe/Amsterdam")
    times2020 = pd.date_range(start='2019-01-01 00:00:00', end='2019-12-31 23:00:00',
                              freq='H', tz=nl_tz)

    # Location of solar panels:
    lon_deg = 5.0
    lat_deg = 52.0
    lon_rad = np.deg2rad(5.0)  # Geographic longitude (>0 for eastern hemisphere; ° -> rad)
    lat_rad = np.deg2rad(52.0)  # Geographic latitude  (>0 for northern hemisphere; ° -> rad)

    # PVLIB
    pos_pv = get_solarposition(times2020, lat_deg, lon_deg, method='nrel_numpy')
    # solarenergy
    times2020_dt = times2020.to_pydatetime()
    sunAz, sunAlt, sunDist = se.sun_position_from_datetime(lon_rad, lat_rad, times2020_dt)
    # sunAz, sunAlt, sunAltUncorr, sunDist = se.sun_position_from_datetime(lon_rad, lat_rad, times2020)

    # calculate differences
    diff_az = (np.rad2deg(sunAz) + 180.0) - pos_pv['azimuth']
    diff_elev = np.rad2deg(sunAlt) - pos_pv['apparent_elevation']

    # plot
    fig, ax = plt.subplots(3, figsize=(15, 8), sharex=True)
    se_style = dict(linestyle='none', marker='o',
                    markerfacecolor='none', markeredgecolor='g', markersize=4)
    pv_style = dict(linestyle='none', marker='.',
                    markerfacecolor='r', markeredgecolor='r', markersize=3)
    ax[0].set_ylabel('Azimuth [$\degree$]')
    ax[0].plot(times2020, pos_pv['azimuth'], label='PVLIB basic', **pv_style)
    ax[0].plot(times2020, np.rad2deg(sunAz) + 180.0, label='SE', **se_style)

    ax[1].set_ylabel('Altitude/apparent elevation [$\degree$]')
    ax[1].plot(times2020, pos_pv['apparent_elevation'], **pv_style)
    ax[1].plot(times2020, np.rad2deg(sunAlt), **se_style)

    # ax[2].plot(times2020, pos1['equation_of_time'])
    ax[2].set_ylabel('Difference [$\degree$]')
    ax[2].plot(times2020, diff_az, label='diff azimuth', **pv_style)
    ax[2].plot(times2020, diff_elev, label='diff alt', **se_style)

    ax[0].legend()
    ax[2].legend()
    ax[2].set_ylim(-1, 1)
    ax[2].set_xlabel('Time')
    plt.suptitle('Validation of SE.sun_position_from_datetime \n'
                 'vs. PVLIB.get_solarposition (apparent_elevation)')
    plt.show()

def test_extra_and_airmass():
    """

    Returns:

    """
    nl_tz = pytz.timezone("Europe/Amsterdam")
    times2020 = pd.date_range(start='2019-01-01 00:00:00', end='2019-12-31 23:00:00',
                              freq='H', tz=nl_tz)

    # Location of solar panels:
    lon_deg = 5.0
    lat_deg = 52.0
    lon_rad = np.deg2rad(5.0)  # Geographic longitude (>0 for eastern hemisphere; ° -> rad)
    lat_rad = np.deg2rad(52.0)  # Geographic latitude  (>0 for northern hemisphere; ° -> rad)

    # PVLIB
    pos_pv = get_solarposition(times2020, lat_deg, lon_deg, method='nrel_numpy')
    # solarenergy
    # times2020_dt = times5060.to_pydatetime()
    sunAz, sunAlt, sunDist = se.sun_position_from_datetime(lon_rad, lat_rad, times2020)

    # PVLIB
    extra_pv = get_extra_radiation(times2020, solar_constant=1361.5, method='asce')
    am_pv = get_relative_airmass(pos_pv['zenith'], model='young1994')

    # solarenergy
    extra_se = se.sol_const / np.square(sunDist)  # Extraterrestrial radiation = Solar constant, scaled with distance
    am_se = se.airmass(sunAlt)                    # Air mass for this Sun altitude

    extFac = se.extinction_factor(am_se)       # Extinction factor at sea level for this air mass
    DNI_se = extra_se / extFac               # DNI for a clear sky

    se_style = dict(linestyle='none', marker='o',
                    markerfacecolor='none', markeredgecolor='g', markersize=4)
    pv_style = dict(linestyle='none', marker='.',
                    markerfacecolor='r', markeredgecolor='r', markersize=3)
    fig, ax = plt.subplots(2,figsize=(15, 8), sharex=True)
    ax[0].set_ylabel('Extra')
    ax[0].plot(times2020, extra_pv, label='PVLIB', **pv_style)
    ax[0].plot(times2020, extra_se, label='SE', **se_style)

    ax[1].set_ylabel('Airmass')
    ax[1].plot(times2020, am_pv, label='PVLIB', **pv_style)
    ax[1].plot(times2020, am_se, label='SE', **se_style)
    plt.suptitle('Validation of Extraterrestrial radiation and Airmass \n '
                 'SE vs. PVLIB')
    plt.show()


def test_distance():
    """

    Returns:

    """
    nl_tz = pytz.timezone("Europe/Amsterdam")
    times2020 = pd.date_range(start='2019-01-01 00:00:00', end='2019-12-31 23:00:00',
                              freq='H', tz=nl_tz)

    # Location of solar panels:
    lon_deg = 5.0
    lat_deg = 52.0
    lon_rad = np.deg2rad(5.0)  # Geographic longitude (>0 for eastern hemisphere; ° -> rad)
    lat_rad = np.deg2rad(52.0)  # Geographic latitude  (>0 for northern hemisphere; ° -> rad)

    # VALIDATION of solar position calculations in PVLIB and solarenergy
    dist_pv = nrel_earthsun_distance(times2020)

    # solarenergy
    #times2020_dt = times5060.to_pydatetime()
    sunAz, sunAlt, sunDist = se.sun_position_from_datetime(lon_rad, lat_rad, times2020)

    # calculate ratio and difference
    ratio_dist = sunDist / dist_pv.values
    diff_dist = sunDist - dist_pv.values

    # index = range(len(sunDist))
    # plot
    fig, ax = plt.subplots(3, figsize=(15, 8), sharex=True)
    se_style = dict(linestyle='none', marker='o',
                    markerfacecolor='none', markeredgecolor='g', markersize=5)
    pv_style = dict(linestyle='none', marker='.',
                    markerfacecolor='r', markeredgecolor='r', markersize=3)
    ax[0].set_ylabel('Distance [A.U.]')
    ax[0].plot(times2020, dist_pv.values, ',r', label='PVLIB basic', **pv_style)
    ax[0].plot(times2020, sunDist, label='SE', **se_style)

    ax[1].set_ylabel('Ratio SE/PVLIB')
    ax[1].plot(times2020, ratio_dist, **pv_style)

    ax[2].set_ylabel('Difference SE - PVLIB')
    ax[2].plot(times2020, diff_dist, **pv_style)

    ax[0].legend()
    ax[2].legend()
    ax[1].set_ylim(0.999, 1.001)
    ax[2].set_ylim(-0.001, 0.001)
    ax[2].set_xlabel('Time')
    plt.suptitle('Validation of AU distance for SE.sun_position_from_datetime \n'
                 'vs. PVLIB.get_solarposition (apparent_elevation)')
    plt.show()


if __name__ == "__main__":
    test_extinction()
    test_positions()
    test_extra_and_airmass()
    test_distance()

