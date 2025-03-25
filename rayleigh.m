n = 100;
[alt, phi] = meshgrid(linspace(0, pi/2, n), linspace(0, 2*pi, 4*n));
theta = pi/2 - alt;

#https://aa.usno.navy.mil/data/AltAz
alt_sun = deg2rad(49.8);
theta_sun = pi/2 - alt_sun;
phi_sun = deg2rad(206.2);

scattering_angle = acos(cos(theta).*cos(theta_sun) + sin(theta_sun)*sin(theta).*cos(phi-phi_sun));
dop = (sin(scattering_angle).^2)./(1 + cos(scattering_angle).^2);

aop = atan((sin(theta)*cos(theta_sun) - cos(theta).*cos(phi-phi_sun)*sin(theta_sun)) ./ sin(phi-phi_sun)*sin(theta_sun));

[x, y, z] = sph2cart(phi, alt, 1);
[x_sun, y_sun, z_sun] = sph2cart(phi_sun, alt_sun, 1);

figure(2);
grid on;
hold on;
scatter3(x_sun, y_sun, z_sun, 200, "filled", "MarkerFaceColor","black")
#surf(x, y, z, aop, "FaceAlpha",0.5, "EdgeAlpha",0)
surf(x, y, z, dop, "FaceAlpha",0.5, "EdgeAlpha",0)
view([-45 45])
colormap jet
colorbar