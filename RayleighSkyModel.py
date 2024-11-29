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

def create_2d_coordinate_set(num_altitude_vals = 91, num_azimuth_vals = 360):
    altitude_observer, azimuth_observer = np.mgrid[0:np.pi/2:num_altitude_vals*1j, 0:np.pi*2:num_azimuth_vals*1j]
    return altitude_observer, azimuth_observer

def create_2d_plot(polar_radius, polar_theta, phi_sun, theta_sun, data, date, cbar_ticks, data_label):
    plt.close("all")
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    #Plot negative azimuth values as the azimuth direction is clockwise but polar plots are counter clockwise
    fig.colorbar(ax.pcolormesh(-polar_theta, polar_radius, data, cmap="jet", vmin=cbar_ticks[0], vmax=cbar_ticks[-1]), ax=ax, ticks=cbar_ticks).set_label(data_label)
    fig.suptitle(date.strftime("%m/%d/%Y %H:%M") + f" UTC | AZ: {round(np.rad2deg(phi_sun),2)}° ALT: {round(np.rad2deg((np.pi/2 - theta_sun)),2)}°")
    ax.grid(False)
    #ax.axis("off")
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.annotate("N", xy=(0,1), xytext=(0,1.08), horizontalalignment="left", verticalalignment="center", arrowprops=dict(facecolor="black",arrowstyle="<|-"))
    ax.annotate("E", xy=(-np.pi/2,1), xytext=(-np.pi/2,1.08), horizontalalignment="center", verticalalignment="top", arrowprops=dict(facecolor="black",arrowstyle="<|-"))
    ax.annotate("S", xy=(-np.pi,1), xytext=(-np.pi,1.08), horizontalalignment="right", verticalalignment="center", arrowprops=dict(facecolor="black",arrowstyle="<|-"))
    ax.annotate("W", xy=(-np.pi*3/2,1), xytext=(-np.pi*3/2,1.08), horizontalalignment="center", verticalalignment="bottom", arrowprops=dict(facecolor="black",arrowstyle="<|-"))
    ax.scatter(-phi_sun, np.tan(np.pi/4 - (np.pi/2 - theta_sun)/2), c="black", s=250, label="Sun", alpha=0.5)
    ax.legend(loc=3, frameon=False)
    fig.tight_layout()
    return fig

def create_multi_plot(phi_sun, theta_sun, date):
    altitude_observer, phi_observer = create_2d_coordinate_set(182, 720)
    theta_observer = np.deg2rad(90) - altitude_observer
    degree_of_polarization = calculate_dop(theta_sun, phi_sun, theta_observer, phi_observer)
    angle_of_polarization = calculate_aop(theta_sun, phi_sun, theta_observer, phi_observer)
    polar_theta = phi_observer
    polar_radius = np.tan(np.pi/4 - altitude_observer/2)
    plt.close("all")
    fig, axis = plt.subplots(1, 2, subplot_kw=dict(projection='polar'))
    dop_ticks = np.linspace(0,100,5)
    aop_ticks = np.linspace(-90,90,5)
    fig.suptitle(f"{date.year}/{date.month}/{date.day} - {date.hour}:{date.minute} UTC | AZ: {round(np.rad2deg(phi_sun),2)} EL: {round(np.rad2deg((np.pi/2 - theta_sun)),2)}")
    fig.colorbar(axis[0].pcolormesh(-polar_theta, polar_radius, degree_of_polarization*100, cmap="jet", vmin=dop_ticks[0], vmax=dop_ticks[-1]), shrink=0.5, ax=axis[0], ticks=dop_ticks, location="left").set_label("Degree of Polarization (%)")
    fig.colorbar(axis[1].pcolormesh(-polar_theta, polar_radius, np.rad2deg(angle_of_polarization), cmap="jet", vmin=aop_ticks[0], vmax=aop_ticks[-1]), shrink=0.5, ax=axis[1], ticks=aop_ticks, location="right").set_label("Angle of Polarization (°)")
    for ax in axis:
        ax.grid(False)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.annotate("N", xy=(0,1), xytext=(0,1.05), horizontalalignment="left", verticalalignment="center")
        ax.annotate("E", xy=(-np.pi/2,1), xytext=(-np.pi/2,1.05), horizontalalignment="center", verticalalignment="top")
        ax.annotate("S", xy=(-np.pi,1), xytext=(-np.pi,1.05), horizontalalignment="right", verticalalignment="center")
        ax.annotate("W", xy=(-np.pi*3/2,1), xytext=(-np.pi*3/2,1.05), horizontalalignment="center", verticalalignment="bottom")
        ax.scatter(-phi_sun, np.tan(np.pi/4 - (np.pi/2 - theta_sun)/2), c="black", s=250, label="Sun", alpha=0.5)
        ax.legend(loc=3, frameon=False)
        fig.tight_layout()
    return fig

def create_2d_plot_dop(phi_sun, theta_sun, date):
    altitude_observer, phi_observer = create_2d_coordinate_set(182, 720)
    theta_observer = np.deg2rad(90) - altitude_observer
    degree_of_polarization = calculate_dop(theta_sun, phi_sun, theta_observer, phi_observer)
    polar_theta = phi_observer
    polar_radius = np.tan(np.pi/4 - altitude_observer/2)
    return create_2d_plot(polar_radius, polar_theta, phi_sun, theta_sun, degree_of_polarization*100, date, np.linspace(0,100,5), "Degree of Polarization (%)")

def create_2d_plot_aop(phi_sun, theta_sun, date):
    altitude_observer, phi_observer = create_2d_coordinate_set(182, 720)
    theta_observer = np.deg2rad(90) - altitude_observer
    angle_of_polarization = calculate_aop(theta_sun, phi_sun, theta_observer, phi_observer)
    polar_theta = phi_observer
    polar_radius = np.tan(np.pi/4 - altitude_observer/2)
    return create_2d_plot(polar_radius, polar_theta, phi_sun, theta_sun, np.rad2deg(angle_of_polarization), date, np.linspace(-90,90,5), "Angle of Polarization (°)")

def create_3d_coordinate_set():
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