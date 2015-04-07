file = 'Latmix2011Site1Profile.nc';
z = ncread(file,'z');
rho = ncread(file,'rho');

latitude = 31;
width = 10e3;
height = 10e3;
Nx = 256;
Nz = 512;
L_z = z(end)-z(1);
z_cheb = abs(L_z/2)*(cos(((0:Nz-1)')*pi/(Nz-1))+1) + min(z);
f0 = 2*(7.2921e-5)*sin(latitude*pi/180);

rho_cheb = interp1(z,rho,z_cheb);
figure
plot(rho,z,'b')
hold on
plot(rho_cheb,z_cheb,'r')

dx = width/Nx;
k_max = 2*pi/dx;

[F0, G0, h0, N2] = InternalWaveModesFromDensityProfile_Spectral( rho_cheb, z_cheb, z_cheb, 0.0, latitude, 'max_u', 'rigid_lid' );

% Note that if you switch this to 'rigid_lid', then the counting of the
% modes will change. The mode in the first position will be the first
% baroclinic mode, rather than the barotropic mode.
[F_max, G_max, h_max, ~] = InternalWaveModesFromDensityProfile_Spectral( rho_cheb, z_cheb, z_cheb, k_max, latitude, 'max_u', 'rigid_lid' );
omega_max = (sqrt(9.81*h_max*k_max*k_max+f0*f0)/(2*pi))';

figure('Position', [0 0 1500 800])
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

subplot(1,2,1)
plot(sqrt(N2)/(2*pi), z_cheb)
ylim([min(z_cheb) max(z_cheb)])
xlabel('buoyancy frequency (cycles per second)')
ylabel('depth (meters)')
hold on
plot(omega_max(1)*ones(size(z_cheb)),z_cheb)


subplot(1,2,2)
plot( F0(:,1), z_cheb, 'LineWidth',2);
ylim([min(z_cheb) max(z_cheb)])
hold on
plot( F_max(:,1), z_cheb, 'LineWidth',2);
%title(sprintf('First two baroclinic (u,v)-modes, F(z), h = %.2g, %.2g, %.2g.', h(1), h(2), h(3) ));
%title(sprintf('First three baroclinic (u,v)-modes, F(z), h = %.2g, %.2g, %.2g.\nBarotropic h = %4g.', h(2), h(3), h(4), h(1) ));
ylabel('Depth');
xlabel('Mode');