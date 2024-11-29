from RayleighSkyModel import *

lat = 38.8976841639134
lon = 77.03653291227445

#lat = 00.0
#lon = 00.0

for h in np.arange(0, 24):
    for m in np.arange(0, 60, 10):
        time = Time(f"2024-9-21 {h}:{m}:00")
        date = time.to_datetime()
        altitude_sun, phi_sun = np.deg2rad(get_solar_coordinates(time, lat, lon))
        theta_sun = np.deg2rad(90) - altitude_sun

        timestamp = f"{h}h{m}m"

        print(timestamp)
        print("Altitude ", np.rad2deg(altitude_sun))
        print("Azimuth ", np.rad2deg(phi_sun))

        if altitude_sun < 0:
            continue

        dop_plot = create_2d_plot_dop(phi_sun, theta_sun, date)
        aop_plot = create_2d_plot_aop(phi_sun, theta_sun, date)

        aop_plot.savefig(f"out\\2d\\aop\\{timestamp}.png", dpi=250)
        dop_plot.savefig(f"out\\2d\\dop\\{timestamp}.png", dpi=250)

create_gif("out\\2d\\dop\\", "2d_dop")
create_gif("out\\2d\\aop\\", "2d_aop")
