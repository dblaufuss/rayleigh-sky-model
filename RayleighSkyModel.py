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

def get_solar_coordinates(time: Time, latitude, longitude):
    sunaltaz = get_sun(time).transform_to(AltAz(obstime=time, location=EarthLocation(lat=latitude*u.deg, lon=longitude*u.deg)))
    return (sunaltaz.alt*u.deg).value, (sunaltaz.az*u.deg).value

def calculate_dop(theta_sun, phi_sun, theta_observer, phi_observer):
    scattering_angle = np.arccos(np.cos(theta_observer)*np.cos(theta_sun) + np.sin(theta_observer)*np.sin(theta_sun)*np.cos(phi_observer - phi_sun))
    return np.power(np.sin(scattering_angle), 2)/(1 + np.power(np.cos(scattering_angle), 2))

def calculate_aop(theta_sun, phi_sun, theta_observer, phi_observer):
    return np.arctan((np.sin(theta_observer)*np.cos(theta_sun) - np.cos(theta_observer)*np.cos(phi_observer - phi_sun)*np.sin(theta_sun)) 
        / np.sin(phi_observer - phi_sun)*np.sin(theta_sun))

def create_2d_coordinate_set(num_altitude_vals = 90, num_azimuth_vals = 360):
    altitude_observer, azimuth_observer = np.mgrid[0:np.pi/2:num_altitude_vals+1j, 0:np.pi*2:num_azimuth_vals*1j]
    return altitude_observer, azimuth_observer

def create_2d_plot(polar_radius, polar_theta, phi_sun, theta_sun, data, date, cbar_ticks):
    plt.close("all")
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    fig.colorbar(ax.pcolormesh(polar_theta, polar_radius, data, cmap="jet", vmin=cbar_ticks[0], vmax=cbar_ticks[-1]), ax=ax, ticks=cbar_ticks)
    ax.set_title(f"{date.year}/{date.month}/{date.day} - {date.hour}:{date.minute} UTC")
    ax.grid(False)
    ax.axis("off")
    ax.scatter(phi_sun, np.tan(np.pi/4 - (np.pi/2 - theta_sun)/2), c="black", s=250, label="Sun", alpha=0.5)
    ax.legend(loc=3, frameon=False)
    return fig

def create_2d_plot_dop(phi_sun, theta_sun, date):
    altitude_observer, phi_observer = create_2d_coordinate_set()
    theta_observer = np.deg2rad(90) - altitude_observer
    degree_of_polarization = calculate_dop(theta_sun, phi_sun, theta_observer, phi_observer)
    polar_theta = phi_observer
    polar_radius = np.tan(np.pi/4 - altitude_observer/2)
    return create_2d_plot(polar_radius, polar_theta, phi_sun, theta_sun, degree_of_polarization, date, np.linspace(0,1,5))

def create_2d_plot_aop(phi_sun, theta_sun, date):
    altitude_observer, phi_observer = create_2d_coordinate_set()
    theta_observer = np.deg2rad(90) - altitude_observer
    angle_of_polarization = calculate_aop(theta_sun, phi_sun, theta_observer, phi_observer)
    polar_theta = phi_observer
    polar_radius = np.tan(np.pi/4 - altitude_observer/2)
    return create_2d_plot(polar_radius, polar_theta, phi_sun, theta_sun, np.rad2deg(angle_of_polarization), date, np.linspace(-90,90,5))

def create_3d_coordinate_set(camera_fov):
    altitude_observer, azimuth_observer = np.mgrid[0:np.pi/2:270j, 0:np.pi*2:720j]
    x_observer = np.cos(altitude_observer)*np.cos(azimuth_observer)
    y_observer = np.cos(altitude_observer)*np.sin(azimuth_observer)
    z_observer = np.sin(altitude_observer)
    return altitude_observer, azimuth_observer, x_observer, y_observer, z_observer

def create_3d_plot(x, y, z, data, date, name):
    fig = go.Figure(go.Surface(x=x, y=y, z=z, surfacecolor=data, colorscale="jet"))
    fig.update_layout(
        title=f"{name} | {date.year}/{date.month}/{date.day} - {date.hour}:{date.minute} UTC",
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
    return fig

def create_3d_plot_dop(phi_sun, theta_sun, date):
    altitude_observer, azimuth_observer, x_observer, y_observer, z_observer = create_3d_coordinate_set()
    phi_observer = azimuth_observer
    theta_observer = np.pi/2 - altitude_observer
    degree_of_polarization = calculate_dop(theta_sun, phi_sun, theta_observer, phi_observer)
    return create_3d_plot(x_observer, y_observer, z_observer, degree_of_polarization, date, "Degree of Polarization")

def create_3d_plot_aop(phi_sun, theta_sun, date):
    altitude_observer, azimuth_observer, x_observer, y_observer, z_observer = create_3d_coordinate_set()
    phi_observer = azimuth_observer
    theta_observer = np.pi/2 - altitude_observer
    angle_of_polarization = calculate_aop(theta_sun, phi_sun, theta_observer, phi_observer)
    return create_3d_plot(x_observer, y_observer, z_observer, angle_of_polarization, date, "Angle of Polarization")

def create_gif(path, name):
    with imageio.get_writer(f"{name}.gif", mode="I") as writer:
        directory = os.listdir(path)
        directory.sort(key=lambda x: os.path.getmtime(path+x))
        for file in directory:
            image = imageio.imread(path + file)
            writer.append_data(image)