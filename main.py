from astropy.coordinates import *
from astropy.time import Time
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib.colors
import numpy as np
import datetime

def gen_rayleigh_model_dop(phi_sun, theta_sun, date, name):
    theta_range = np.linspace(0, 2*np.pi, 1080)
    r_range = np.linspace(0, 1, 1000)
    polar_radius, polar_theta = np.meshgrid(r_range, theta_range)

    phi_observer = polar_theta
    theta_observer = np.arcsin(polar_radius)

    scattering_angle = np.arccos(np.cos(theta_observer)*np.cos(theta_sun) + np.sin(theta_observer)*np.sin(theta_sun)*np.cos(phi_observer - phi_sun))
    degree_of_polarization = np.power(np.sin(scattering_angle), 2)/(1 + np.power(np.cos(scattering_angle), 2))

    plt.subplot(projection="polar")
    plt.suptitle(f"{date.year}/{date.month}/{date.day} - {date.hour}:{date.minute} UTC")
    plt.pcolormesh(polar_theta, polar_radius, degree_of_polarization, cmap="jet")
    plt.colorbar()
    plt.grid(False)
    plt.axis("off")
    plt.scatter(phi_sun, np.sin(theta_sun), c="black", s=250, label="Sun", alpha=0.5)
    plt.legend(loc=3, frameon=False)
    plt.savefig(f"out\\dop\\{name}.png", dpi=500)


def gen_rayleigh_model_aop(phi_sun, theta_sun, date, name):
    theta_range = np.linspace(0, 2*np.pi, 1080)
    r_range = np.linspace(0, 1, 1000)
    polar_radius, polar_theta = np.meshgrid(r_range, theta_range)

    phi_observer = polar_theta
    theta_observer = np.arcsin(polar_radius)

    angle_of_polarization = np.arctan(
        np.sin(theta_observer)*np.cos(theta_sun) - np.cos(theta_observer)*np.cos(phi_observer - phi_sun)*np.sin(theta_sun) /
        np.sin(phi_observer - phi_sun)*np.sin(theta_sun)
    )

    plt.subplot(projection="polar")
    plt.suptitle(f"{date.year}/{date.month}/{date.day} - {date.hour}:{date.minute} UTC")
    plt.pcolormesh(polar_theta, polar_radius, np.rad2deg(angle_of_polarization), cmap="jet")
    plt.colorbar()
    plt.grid(False)
    plt.axis("off")
    plt.scatter(phi_sun, np.sin(theta_sun), c="black", s=250, label="Sun", alpha=0.5)
    plt.legend(loc=3, frameon=False)
    plt.savefig(f"out\\aop\\{name}.png", dpi=500)

for h in np.arange(0, 24):
    for m in np.arange(0, 60, 10):
        time = Time(f"2024-10-23 {h}:{m}:00")
        date = time.to_datetime(timezone=datetime.timezone.utc)
        sun = get_sun(time)
        observer = EarthLocation(lat=38.8976841639134*u.deg, lon=-77.03653291227445*u.deg)
        sunaltaz = sun.transform_to(AltAz(obstime=time, location=observer))

        phi_sun = np.deg2rad((sunaltaz.az*u.deg).value)
        altitude_sun = np.deg2rad((sunaltaz.alt*u.deg).value)
        theta_sun = np.deg2rad(90) - altitude_sun

        name = f"{h}h{m}m"

        print(name)
        print("Azimuth ", (sunaltaz.az*u.deg).value)
        print("Altitude ", (sunaltaz.alt*u.deg).value)

        if altitude_sun < 0:
            continue

        gen_rayleigh_model_dop(phi_sun, theta_sun, date, name)
        gen_rayleigh_model_aop(phi_sun, theta_sun, date, name)