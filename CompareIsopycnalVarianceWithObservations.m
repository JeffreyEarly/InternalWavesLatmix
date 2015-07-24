file = '/Volumes/home/jearly/InternalWavesLatmix_256_256_80_GM_0.031.nc';
observations = open('Site1GliderIsopycnalVariance.mat');

x = ncread(file, 'x');
y = ncread(file, 'y');
z = ncread(file, 'z');
N2 = ncread(file, 'N2');

% Extract and reshape zeta so that it's depth x samples
xStride = 8;
zeta = double(squeeze(ncread(file, 'zeta', [1 1 1 1], [length(y)/xStride length(x)/xStride length(z) 1], [xStride xStride 1 1])));
zeta = (reshape(zeta,[size(zeta,1)*size(zeta,2) size(zeta,3)]))';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Plot the WKB scale isopycnal variance as a function of depth.
%
N0 = 5.2E-3; % radians per second, GM frequency
b = 1300; % meters, GM depth
figure
plot(squeeze(mean(zeta.^2,2)).*sqrt(N2)/N0,z,'k','LineWidth',2)
hold on
plot(observations.wkbIsopycnalVariance, observations.depth,'blue','LineWidth',2);
vlines(53, 'g--')
ylim([-100 0])
legend('Model', 'Glider Observations', 'GM')
xlabel('Isopycnal displacement variance (m^2)')
ylabel('depth (m)')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Plot spectrum of the isopycnal displacement
%

% First compute the GM version of the spectrum
j_star = 3*pi/b;
H = (3/2)*(j_star^(3/2))*(kz + j_star).^(-5/2);
gm = (1/2)*(b*b*6.3E-5)*H;

% Project the model data onto a uniform grid
zUniform = (z(1):(min(diff(z))):z(end))';
zetaUniform = interp1(z,zeta,zUniform);

dz = abs(zUniform(2)-zUniform(1));
taper_bandwidth = 1;
[psi,lambda]=sleptap(size(zetaUniform,1),taper_bandwidth);
[kz,s]=mspec(dz,zetaUniform,psi);
s = s/(2*pi); % force it to sum to variance

% Limit the range to the largest wavenumber of the original non-uniform
% grid
validIndices = find(kz/(2*pi) < 1/(2*max(diff(z))));
kz = kz(validIndices);
s = s(validIndices,:);

figure
% plot(kz/(2*pi),s)
% hold on,
plot(kz/(2*pi),mean(s,2),'k','LineWidth',2)
hold on,
plot(observations.kz/(2*pi), observations.Smean,'blue','LineWidth',2);
plot(kz/(2*pi),gm,'g--')
ylog, xlog
legend('Model', 'Glider Observations', 'GM')
xlabel('frequency (cycles per meter)')
ylabel('power (m^3)')
xlim([kz(2)/(2*pi) max(kz/(2*pi))])

% y = log10(mean(s,2)); y=y(2:end);
% x = log10(kz); x=x(2:end);
% [p,S,mu]=polyfit(x,y,1);
% slope_fit = p(1)/mu(2);
% intercept_fit = p(2) - p(1)*mu(1)/mu(2);
% s_fit = (10^intercept_fit)*kz.^slope_fit;
% plot(kz/(2*pi),s_fit,'g','LineWidth',2)

