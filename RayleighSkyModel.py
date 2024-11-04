from astropy.coordinates import *
from astropy.time import Time
from datetime import datetime
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import imageio.v2 as imageio
import os

def get_solar_position(time: Time, latitude, longitude):
    sun = get_sun(time)
    observer = EarthLocation(lat=latitude*u.deg, lon=longitude*u.deg)
    sunaltaz = sun.transform_to(AltAz(obstime=time, location=observer))
    altitude = (sunaltaz.alt*u.deg).value
    azimuth = (sunaltaz.az*u.deg).value
    return altitude, azimuth

def calculate_dop(theta_sun, phi_sun, theta_observer, phi_observer):
    scattering_angle = np.arccos(np.cos(theta_observer)*np.cos(theta_sun) + np.sin(theta_observer)*np.sin(theta_sun)*np.cos(phi_observer - phi_sun))
    degree_of_polarization = np.power(np.sin(scattering_angle), 2)/(1 + np.power(np.cos(scattering_angle), 2))
    return degree_of_polarization

def calculate_aop(theta_sun, phi_sun, theta_observer, phi_observer):
    angle_of_polarization = np.arctan(
        (np.sin(theta_observer)*np.cos(theta_sun) - np.cos(theta_observer)*np.cos(phi_observer - phi_sun)*np.sin(theta_sun))
        / np.sin(phi_observer - phi_sun)*np.sin(theta_sun)
    )
    return angle_of_polarization

def create_2d_coordinate_set(num_theta_vals, num_radii_vals):
    theta_range = np.delete(np.linspace(0, 2*np.pi, num_theta_vals), -1)
    radius_range = np.linspace(0, 1, 500)
    polar_radius, polar_theta = np.meshgrid(radius_range, theta_range)
    return polar_radius, polar_theta

def create_2d_plot_dop(phi_sun, theta_sun, date, name):
    polar_radius, polar_theta = create_2d_coordinate_set(360, 500)

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
    plt.savefig(f"out\\2d\\dop\\{name}.png", dpi=250)

def create_2d_plot_aop(phi_sun, theta_sun, date, name):
    polar_radius, polar_theta = create_2d_coordinate_set(360, 500)
    phi_observer = polar_theta
    theta_observer = np.arcsin(polar_radius)
    angle_of_polarization = calculate_aop(theta_sun, phi_sun, theta_observer, phi_observer)

    plt.subplot(projection="polar")
    plt.suptitle(f"{date.year}/{date.month}/{date.day} - {date.hour}:{date.minute} UTC")
    plt.pcolormesh(polar_theta, polar_radius, np.rad2deg(angle_of_polarization), cmap="jet")
    plt.colorbar()
    plt.grid(False)
    plt.axis("off")
    plt.scatter(phi_sun, np.sin(theta_sun), c="black", s=250, label="Sun", alpha=0.5)
    plt.legend(loc=3, frameon=False)
    plt.savefig(f"out\\2d\\aop\\{name}.png", dpi=250)

def create_3d_coordinate_set():
    altitude_observer, azimuth_observer = np.mgrid[0:np.pi/2:270j, 0:np.pi*2:720j]

    x_observer = np.cos(altitude_observer)*np.cos(azimuth_observer)
    y_observer = np.cos(altitude_observer)*np.sin(azimuth_observer)
    z_observer = np.sin(altitude_observer)

    return altitude_observer, azimuth_observer, x_observer, y_observer, z_observer

def create_3d_plot_dop(phi_sun, theta_sun, date, name):
    altitude_observer, azimuth_observer, x_observer, y_observer, z_observer = create_3d_coordinate_set()

    phi_observer = azimuth_observer
    theta_observer = np.pi/2 - altitude_observer

    degree_of_polarization = calculate_dop(theta_sun, phi_sun, theta_observer, phi_observer)

    fig = go.Figure(go.Surface(x=x_observer, y=y_observer, z=z_observer, surfacecolor=degree_of_polarization, colorscale="jet"))
    fig.update_layout(
        title=f"Degree of Polarization | {date.year}/{date.month}/{date.day} - {date.hour}:{date.minute} UTC",
        template='simple_white',
        scene = dict(
            xaxis = dict(visible=False),
            yaxis = dict(visible=False),
            zaxis = dict(visible=False)
        ),
                scene_camera = dict(
            up=dict(x=0, y=0, z=1),
            center=dict(x=0, y=0, z=0),
            eye=dict(x=1, y=1, z=1)
        )
    )
    fig.write_image(f"out\\3d\\dop\\{name}.png", engine="kaleido")

def create_3d_plot_aop(phi_sun, theta_sun, date, name):
    altitude_observer, azimuth_observer, x_observer, y_observer, z_observer = create_3d_coordinate_set()

    phi_observer = azimuth_observer
    theta_observer = np.pi/2 - altitude_observer
    
    angle_of_polarization = calculate_aop(theta_sun, phi_sun, theta_observer, phi_observer)

    fig = go.Figure(go.Surface(x=x_observer, y=y_observer, z=z_observer, surfacecolor=np.rad2deg(angle_of_polarization), colorscale="jet"))
    fig.update_layout(
        title=f"Degree of Polarization | {date.year}/{date.month}/{date.day} - {date.hour}:{date.minute} UTC",
        template='simple_white',
        scene = dict(
            xaxis = dict(visible=False),
            yaxis = dict(visible=False),
            zaxis = dict(visible=False)
        ),
                scene_camera = dict(
            up=dict(x=0, y=0, z=1),
            center=dict(x=0, y=0, z=0),
            eye=dict(x=1, y=1, z=1)
        )
    )
    fig.write_image(f"out\\3d\\aop\\{name}.png", engine="kaleido")

def create_gif(path, name):
    with imageio.get_writer(f"{name}.gif", mode="I") as writer:
        directory = os.listdir(path)
        directory.sort(key=lambda x: os.path.getmtime(path+x))

        for file in directory:
            image = imageio.imread(path + file)
            writer.append_data(image)