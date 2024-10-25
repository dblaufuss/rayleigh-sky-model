from pysolar.solar import *
import datetime
import numpy as np
import matplotlib.pyplot as plt

def gen_dop_image(date: datetime.datetime, lat, lon, name):
    phi_sun = np.deg2rad(get_azimuth(lat, lon, date))
    altitude_sun = np.deg2rad(get_altitude(lat, lon, date))
    theta_sun = np.deg2rad(90) - altitude_sun

    theta_vals = []
    r_vals = []
    aop_vals = []
    dop_vals = []

    for theta_observer in np.arange(0, 90.5, 0.5):
        for phi_observer in np.arange(0, 360.5, 0.5):
            theta_vals.append(np.deg2rad(phi_observer))
            r_vals.append(np.sin(np.deg2rad(theta_observer)))

            angle_of_polarization = np.rad2deg(np.arctan2(
                np.sin(np.deg2rad(theta_observer))*np.cos(theta_sun) - np.cos(np.deg2rad(theta_observer))*np.cos(np.deg2rad(phi_observer) - phi_sun)*np.sin(theta_sun),
                np.sin(np.deg2rad(phi_observer) - phi_sun)*np.sin(theta_sun)
            ))
            if angle_of_polarization < -90:
                angle_of_polarization += 180
            elif angle_of_polarization > 90:
                angle_of_polarization -= 180
            aop_vals.append(angle_of_polarization)

            scattering_angle = np.arccos(np.cos(np.deg2rad(theta_observer))*np.cos(theta_sun) + np.sin(np.deg2rad(theta_observer))*np.sin(theta_sun)*np.cos(np.deg2rad(phi_observer) - phi_sun))
            degree_of_polarization = np.power(np.sin(scattering_angle), 2)/(1 + np.power(np.cos(scattering_angle), 2))
            dop_vals.append(degree_of_polarization)

    fig = plt.subplot(projection = "polar")
    plt.suptitle(f"{date.year}/{date.month}/{date.day} - {date.hour}:{date.minute} UTC")
    fig.set_title("Degree of Polarization")
    fig.scatter(theta_vals, r_vals, c = dop_vals, cmap = "jet")
    fig.grid(False)
    fig.axis("off")
    fig.set_xticks([])
    fig.set_yticks([])

    plt.savefig(f"out\\dop\\{name}.png")
    plt.clf()

    fig = plt.subplot(projection = "polar")
    plt.suptitle(f"{date.year}/{date.month}/{date.day} - {date.hour}:{date.minute} UTC")
    fig.set_title("Angle of Polarization")
    fig.scatter(theta_vals, r_vals, c = aop_vals, cmap = "jet")
    fig.grid(False)
    fig.axis("off")
    fig.set_xticks([])
    fig.set_yticks([])

    plt.savefig(f"out\\aop\\{name}.png")
    plt.clf()

latitude = 38.8976841639134
longitude = -77.03653291227445

for h in range(0,23):
    day = datetime.datetime(
        year = 2024,
        month = 10,
        day = 24,
        hour = h,
        minute = 0,
        second = 0,
        tzinfo = datetime.timezone.utc
    )
    gen_dop_image(day, latitude, longitude, h)

day.day = 25
day.hour = 0
gen_dop_image(day, latitude, longitude, 24)


with imageio.get_writer("out\\aop.gif", mode="I") as writer:
    for file in os.listdir("out\\aop\\"):
        image = imageio.imread(file)
        writer.append_data(image)

with imageio.get_writer("out\\dop.gif", mode="I") as writer:
    for file in os.listdir("out\\dop\\"):
        image = imageio.imread(file)
        writer.append_data(image)