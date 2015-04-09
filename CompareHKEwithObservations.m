N0 = 5.2e-3; % reference buoyancy frequency, radians/seconds
N_50 = 2*pi*1.7e-3; % radians per second

% Model
file = '/Volumes/jearly/Desktop/InternalWavesLatmix_256_256_50_GM_0.042.nc';
file = '/Volumes/home/jearly/InternalWavesLatmix_256_256_50_GM_0.062.nc';
% file = '/Volumes/Data/InternalWavesLatmix_256_256_50_GM_0.062.nc';
% file = '/Volumes/Data/InternalWaveSimulations/InternalWavesGMSpectrumExponentialStratification.nc';
N2 = double(ncread(file, 'N2'));
f0 = ncreadatt(file, '/', 'f0');
latitude = ncreadatt(file, '/', 'latitude');
x = ncread(file, 'x');
y = ncread(file, 'y');
z = ncread(file, 'z');
t = ncread(file, 'time');
xStride = 4;

depth = -50;
[depth_index] = find(z >= depth, 1, 'first');


maxTime = 34*3600;
t_index = find(t <= maxTime, 1, 'last');

u3d = double(squeeze(ncread(file, 'u', [1 1 depth_index 1], [length(y)/xStride length(x)/xStride 1 t_index], [xStride xStride 1 1])));
v3d = double(squeeze(ncread(file, 'v', [1 1 depth_index 1], [length(y)/xStride length(x)/xStride 1 t_index], [xStride xStride 1 1])));

dt = t(2)-t(1);

[M, N, K] = size(u3d);

% Compute a few 'mooring' time series
cv_mooring = zeros([K 1]);
subsample = 1;
iMooring = 0;
for i=1:subsample:M
	for j=1:subsample:N
		iMooring = iMooring+1;
		cv_mooring(:,iMooring) = squeeze(u3d(i,j,:) + sqrt(-1)*v3d(i,j,:));
	end
end

taper_bandwidth = 2;
psi=[];
[psi,lambda]=sleptap(size(cv_mooring,1),taper_bandwidth);
[omega_p, Spp, Snn, Spn] = mspec(dt,cv_mooring,psi);

omega = [ -flipud(omega_p(2:end)); omega_p];
S = [flipud(vmean(Snn,2)); vmean(Spp(2:end,:),2)]/(4*pi);
% A factor of 1/(2*pi) so this sums to the variance, then another factor of
% 1/2 so that it's a kinetic energy.

figure

% The N0/N_50 factor WKB normalizes it, and the 2*pi/3600 converts to RC's
% units.
plot( 3600*omega/(2*pi), (N0/N_50)*2*pi*S/3600, 'red', 'LineWidth', 2)


% Observations
load('LatMix_iwsp.mat')
LatMixSite = 1;

hold on
plot(frequency_cph, CCWsp_m2s2_cph_wkb(:,LatMixSite), 'Color', 'green', 'LineWidth', 2)
plot(-frequency_cph, CWsp_m2s2_cph_wkb(:,LatMixSite), 'Color', 'green', 'LineWidth', 2)
plot(gmfreqcph, gmccwcph, 'Color', 'blue', 'LineWidth', 2)
plot(-gmfreqcph, gmcwcph, 'Color', 'blue', 'LineWidth', 2)
ylog
xlim([-0.5 0.5])
ylim([1e-4 .8e-1])


freq = [-fliplr(gmfreqcph) gmfreqcph];

[S_gm] = GarrettMunkHorizontalKineticEnergyRotarySpectrumWKB( 2*pi*freq/3600, latitude, N0, 0 );
%[S_gm2] = GarrettMunkHorizontalKineticEnergyRotarySpectrum( 2*pi*freq/3600, latitude, N_50, 0 );


plot(freq, 2*pi*S_gm/3600, 'Color', 'black', 'LineWidth', 2)
ylog

legend('Observations', 'RC GM', 'JJE GM WKB')


% Extract the GM Spectrum at the reference level
[S_gm] = GarrettMunkHorizontalKineticEnergySpectrumWKB( 2*pi*gmfreqcph/3600, latitude, N0, 0 );

% Grab the energy at 
[psi,lambda]=sleptap(size(cv_mooring,1),taper_bandwidth);
[omega, Su] = mspec(dt,real(cv_mooring),psi);
[omega, Sv] = mspec(dt,imag(cv_mooring),psi);
HKE_model = vmean((Su + Sv)/(4*pi),2); % A factor of 1/2 to go from variance to KE, and a factor of 2*pi to deal with summation

figure
plot( frequency_cph, HKEsp_m2s2_cph_wkb(:,LatMixSite), 'Color', 'green', 'LineWidth', 2), ylog
hold on
plot( gmfreqcph, gmHKEcph, 'Color', 'blue', 'LineWidth', 2)
plot( gmfreqcph, 2*pi*S_gm/3600, 'Color', 'black', 'LineWidth', 2)
plot( 3600*omega/(2*pi), (N0/N_50)*2*pi*HKE_model/3600, 'Color', 'magenta', 'LineWidth', 2)
xlim([0.0 0.6])
ylim([1e-4 2e-2])