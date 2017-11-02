%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Actual measured profile
%

file = '../Latmix2011Site1Profile.nc';
rhoLatmix = ncread(file,'rho');
N2Latmix = ncread(file,'N2');
zLatmix = ncread(file,'z');

z = linspace(min(zLatmix),max(zLatmix),250);
z = linspace(-300,0,1024);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Realistic stratification structure function
%
rho0 = min(rhoLatmix);
g = 9.81;
latitude = 31;

delta_p = 0.9;
z_s = -100;
L_s = 10;
L_d = 150; %L_d = 1300;
z_p = -17; %z_p = -1;
D = 5000;

N0 = 2.8e-3; % Surface
Nq = 1.4e-2; % Pycnocline portion to exponentials
Np = 4.7e-2;
Nd = 1.75e-3;

if 0
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
elseif 0
    % Notes for this profile are here:
    % /Users/jearly/Dropbox/Documents/Notes/latmix-glider-profiles/ADifferentAnalyticalProfile/ADifferentAnalyticalProfile.tex
    L_s = 8;
    L_d = 100;
    
    A = Np*Np - Nq*Nq;
    B = (Nq*Nq - N0*N0)/(1 - exp(-2*z_p^2/L_s^2));
    C = N0*N0-B*exp(-2*z_p^2/L_s);
    E = (Nq*Nq - Nd*Nd)/( 1 - exp(-2*(-D-z_p)^2/L_d^2) );
    F = Nd*Nd - E*exp( -2*(-D-z_p)^2/L_d^2 );
    
    rho_surface = @(z) rho0*(1 - (L_s*B/(2*g)) * sqrt(pi/2) *( erf(sqrt(2)*z_p/L_s) - erf(-sqrt(2)*(z-z_p)/L_s) ) - C*z/g);
    rho_deep = @(z) rho_surface(z_p) - (rho0*L_d*E/(2*g)) * sqrt(pi/2) * (erf(sqrt(2)*(z-z_p)/L_d)) - rho0*F*(z-z_p)/g;
    rho_p = @(z) -(A*rho0*delta_p/g)*(tanh( (z-z_p)/delta_p) - 1);
    
    rho = @(z) (z>=z_p).*rho_surface(max(z,z_p)) + (z<z_p).*rho_deep(z) + rho_p(z);
    
    N2_surface = @(z) B*exp(-2*(z-z_p).^2/L_s^2) + C;
    N2_deep = @(z) E*exp(-2*(z-z_p).^2/L_d^2) + F;
    N2_p = @(z) A*sech( (z-z_p)/delta_p ).^2;
    
    N2 = @(z) (z>=z_p).*N2_surface(z) + (z<z_p).*N2_deep(z) + N2_p(z);
else
    % Notes for this profile are here:
    % /Users/jearly/Dropbox/Documents/Notes/latmix-glider-profiles/ADifferentAnalyticalProfile/ADifferentAnalyticalProfile.tex
    delta_p = 0.9;
    L_s = 8;
    L_d = 100;
    z_p = -17;
    z_T = -190;
    D = 5000;
    b = 1300;
    
    N0 = 2.8e-3; % Surface
    Nq = 1.4e-2; % Pycnocline portion to exponentials
    Np = 4.7e-2;
    
    A = Np*Np - Nq*Nq;
    B = (Nq*Nq - N0*N0)/(1 - exp(-2*z_p^2/L_s^2));
    C = N0*N0-B*exp(-2*z_p^2/L_s);
    alpha = -(2*b*(z_T-z_p)/(L_d^2))*exp(-2*(z_T-z_p)^2/(L_d^2) - 2*z_T/b);
    gamma = alpha*exp(2*z_T/b)-exp(-2*(z_T-z_p)^2/(L_d^2));
    E = Nq*Nq/( 1 + gamma );
    F = E*gamma;
    G = E*alpha;
     
    rho_surface = @(z) rho0*(1 - (L_s*B/(2*g)) * sqrt(pi/2) *( erf(sqrt(2)*z_p/L_s) - erf(-sqrt(2)*(z-z_p)/L_s) ) - C*z/g);
    rho_mid = @(z) rho_surface(z_p) - (rho0*L_d*E/(2*g)) * sqrt(pi/2) * (erf(sqrt(2)*(z-z_p)/L_d)) - rho0*F*(z-z_p)/g;
    rho_deep = @(z) rho_mid(z_T) - (rho0*b/(2*g))*G*(exp(2*z/b)-exp(2*z_T/b));
    rho_p = @(z) -(A*rho0*delta_p/g)*(tanh( (z-z_p)/delta_p) - 1);
    
    rho = @(z) (z>=z_p).*rho_surface(max(z,z_p)) + (z<z_p & z > z_T).*rho_mid(z) + (z<=z_T).*rho_deep(z) + rho_p(z);
    
    N2_surface = @(z) B*exp(-2*(z-z_p).^2/L_s^2) + C;
    N2_mid = @(z) E*exp(-2*(z-z_p).^2/L_d^2) + F;
    N2_deep = @(z) G*exp(2*z/b);
    N2_p = @(z) A*sech( (z-z_p)/delta_p ).^2;
    
    N2 = @(z) (z>=z_p).*N2_surface(z) + (z<z_p & z >z_T).*N2_mid(z) + (z<=z_T).*N2_deep(z) + N2_p(z);
end

im = InternalModes(rho,[-D 0],z,latitude,  'method', 'wkbSpectral', 'nEVP', 513); % adding extra points to just confirm all is well converged.

fprintf('%f times GM at the bottom.\n',sqrt(N2(-D))/(5.2e-3 * exp(-5000/1300)));

figure
subplot(1,2,1)
plot(N2Latmix,zLatmix),xlog, hold on
plot(N2(z),z)
plot(im.N2,z)
legend('latmix','analytic', 'computed')

subplot(1,2,2)
plot(rhoLatmix,zLatmix), hold on
plot(rho(z),z)

im.ShowLowestModesAtWavenumber(0.0)

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Now let's create a stretched coordinate and see how well it does
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nz = 64;
maxDepth = -50;
z = linspace(0,maxDepth,1e5);

N_scaled = sqrt(N2(z)/N2(0));
s = cumtrapz( z, N_scaled );
sGrid = linspace(min(s), max(s), Nz);
z_stretched = interp1(s, z, sGrid );

nGrid = 2^14+1; % 2^13 looks fine, but 2^14 looks perfect.

method = 'densitySpectral';
im1024 = InternalModes(rho,[-D 0],z_stretched,latitude,  'method', method, 'nEVP', 1024, 'nGrid', nGrid);
im512 = InternalModes(rho,[-D 0],z_stretched,latitude,  'method', method, 'nEVP', 512, 'nGrid', nGrid);
im256 = InternalModes(rho,[-D 0],z_stretched,latitude,  'method', method, 'nEVP', 256, 'nGrid', nGrid);
im128 = InternalModes(rho,[-D 0],z_stretched,latitude,  'method', method, 'nEVP', 128, 'nGrid', nGrid);

% Biggest issue occurs for long wavelength (k=0)...
k = 0;

[F1024,G1024,h1024] = im1024.ModesAtWavenumber(k);
[F512,G512,h512] = im512.ModesAtWavenumber(k);
[F256,G256,h256] = im256.ModesAtWavenumber(k);
[F128,G128,h128] = im128.ModesAtWavenumber(k);

%...and the lowest mode
iMode = 1;
figure
plot(F1024(:,iMode),z_stretched), hold on
plot(F512(:,iMode),z_stretched)
plot(F256(:,iMode),z_stretched)
plot(F128(:,iMode),z_stretched)

legend('1024','512','256','128')

% Looks like you can get away with 512, although 1024 would be better.

