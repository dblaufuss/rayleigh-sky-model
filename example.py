from RayleighSkyModel import *

lat = 38.8976841639134
lon = 77.03653291227445

for h in np.arange(0, 24):
    for m in np.arange(0, 60, 10):
        time = Time(f"2024-10-23 {h}:{m}:00")
        date = time.to_datetime()
        altitude_sun, phi_sun = np.deg2rad(get_solar_position(time, lat, lon))
        theta_sun = np.deg2rad(90) - altitude_sun

        name = f"{h}h{m}m"

        print(name)
        print("Altitude ", np.rad2deg(altitude_sun))
        print("Azimuth ", np.rad2deg(phi_sun))

        if altitude_sun < 0:
            continue

        create_3d_plot_aop(phi_sun, theta_sun, date, name)

create_gif("out\\3d\\aop\\", "3d_aop")