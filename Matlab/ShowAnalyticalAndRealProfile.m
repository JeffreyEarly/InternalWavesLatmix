%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Actual measured profile
%

file = '../Latmix2011Site1Profile.nc';
rhoLatmix = ncread(file,'rho');
N2Latmix = ncread(file,'N2');
zLatmix = ncread(file,'z');

z = linspace(min(zLatmix),max(zLatmix),250);
z = linspace(-1000,0,1024);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Realistic stratification structure function
%
rho0 = 1025;
g = 9.81;
latitude = 31;

delta_p = 2;
z_s = -100;
L_s = 10;
L_d = 150; %L_d = 1300;
z_p = -17; %z_p = -1;
D = 1000;

N0 = 2.8e-3; % Surface
Nq = 1.2e-2; % Pycnocline portion to exponentials
Np = 4.7e-2;
Nd = 2e-3;

if 1
    A = Np*Np - Nq*Nq;
    B = (Nq*Nq- N0*N0)/(exp(-2*(z_p-z_s)/L_s) - exp(2*z_s/L_s));
    C = N0*N0-B*exp(2*z_s/L_s);
    E = (Nd*Nd - Nq*Nq)/( exp(2*(-D-z_p)/L_d) - 1);
    F = Nq*Nq - E;
    
    rho_surface = @(z) rho0*(1 - (L_s*B/(2*g)) * (exp(2*z_s/L_s) - exp(-2*(z-z_s)/L_s)) + C*(z_p-z)/g);
    rho_deep = @(z) rho_surface(z_p) + (rho0*L_d*E/(2*g))*(1-exp(2*(z-z_p)/L_d)) - rho0*F*z/g;
    rho_p = @(z) -(A*rho0*delta_p/g)*(tanh( (z-z_p)/delta_p) - 1);
    
    rho = @(z) (z>=z_p).*rho_surface(max(z,z_p)) + (z<z_p).*rho_deep(z) + rho_p(z);
    
    N2_surface = @(z) B*exp(-2*(z-z_s)/L_s) + C;
    N2_deep = @(z) E*exp(2*(z-z_p)/L_d) + F;
    N2_p = @(z) A*sech( (z-z_p)/delta_p ).^2;
    
    N2 = @(z) (z>=z_p).*N2_surface(z) + (z<z_p).*N2_deep(z) + N2_p(z);
else
    
    L_d = 150;
    
    A = Np*Np - Nq*Nq;
    B = (Nq*Nq - N0*N0)/(1 - exp(-2*z_p^2/L_s^2));
    C = N0*N0-B*exp(-2*z_p^2/L_s);
    E = (Nq*Nq - Nd*Nd)/( 1 - exp(-2*(-D-z_p)^2/L_d^2) );
    F = Nd*Nd - E*exp( -2*(-D-z_p)^2/L_d^2 );
    
%     rho_surface = @(z) rho0*(1 - (L_s*B/(2*g)) * (exp(2*z_s/L_s) - exp(-2*(z-z_s)/L_s)) - C*z/g);
%     rho_deep = @(z) rho_surface(z_p) + (rho0*L_d*E/(2*g))*(1-exp(2*(z-z_p)/L_d)) - rho0*F*z/g;
%     rho_p = @(z) -(A*rho0*delta_p/g)*(tanh( (z-z_p)/delta_p) - 1);
%     
%     rho = @(z) (z>=z_p).*rho_surface(max(z,z_p)) + (z<z_p).*rho_deep(z) + rho_p(z);
    
    N2_surface = @(z) B*exp(-2*(z-z_p).^2/L_s^2) + C;
    N2_deep = @(z) E*exp(-2*(z-z_p).^2/L_d^2) + F;
    N2_p = @(z) A*sech( (z-z_p)/delta_p ).^2;
    
    N2 = @(z) (z>=z_p).*N2_surface(z) + (z<z_p).*N2_deep(z) + N2_p(z);
end


figure
subplot(1,2,1)
plot(N2Latmix,zLatmix), hold on
plot(N2(z),z)
xlog



subplot(1,2,2)
plot(rho(z),z)


im = InternalModes(rho,[-D 0],z,latitude,  'method', 'spectral'); % adding extra points to just confirm all is well converged.

% 'nEVP', 1024,
figure
plot(N2Latmix,zLatmix), hold on
plot(im.N2,z)
xlog

im.ShowLowestModesAtWavenumber(0.0)
