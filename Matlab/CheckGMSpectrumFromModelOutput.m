%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CheckGMSpectrumFromModelOutput
%
% This scripts reads the model output generated from
% SaveModelOutputToNetCDF.m and does a few basics tests to see if the wave
% field is as expected.
%
% Jeffrey J. Early
% jeffrey@jeffreyearly.com
%
% December 6th, 2016      Version 1.0

file = '/Volumes/OceanTransfer/FloatsWithTemperatureProfileExperiment_2017-07-26T210656_256x64x16.nc';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Read in the problem dimensions and parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = ncread(file, 'x');
y = ncread(file, 'y');
z = ncread(file, 'z');
t = ncread(file, 't');

latitude = ncreadatt(file, '/', 'latitude');
N2 = ncread(file, 'N2');
N = sqrt(N2);
N0_GM = 5.2e-3;
f0 = 2*(7.2921e-5)*sin(latitude*pi/180);
Nmax = max(sqrt(N2));

shouldWKBScale = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Check the frequency spectrum at a given depth, but multiple (x,y)
% locations.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

depth = -30;
[depth_index] = find(z <= depth, 1, 'last');

% This is what is wonderful about NetCDF. We can pull out a slice in time
% at a given depth, for a range of (strided!) x,y values.
stride = 4;
t_index = length(t)-1;
u3d = double(squeeze(ncread(file, 'u', [1 1 depth_index 1], [length(x)/stride length(y)/stride 1 t_index], [stride stride 1 1])));
v3d = double(squeeze(ncread(file, 'v', [1 1 depth_index 1], [length(x)/stride length(y)/stride 1 t_index], [stride stride 1 1])));

% Create a few 'mooring' time series from this.
[M, Q, K] = size(u3d);
cv_mooring = zeros([K 1]);
subsample = 1;
iMooring = 0;
for i=1:subsample:M
	for j=1:subsample:Q
		iMooring = iMooring+1;
		cv_mooring(:,iMooring) = squeeze(u3d(i,j,:) + sqrt(-1)*v3d(i,j,:));
	end
end

% Compute the spectrum
[omega_p, Spp, Snn, Spn] = mspec(t(2)-t(1),cv_mooring,[]);
omega = [ -flipud(omega_p(2:end)); omega_p];
% We want the integral of this to give us the variance back, so we need to
% divide by 2*pi
S = (1/(2*pi))*[flipud(vmean(Snn,2)); vmean(Spp(2:end,:),2)];

figure, plot(omega,S), ylog

xAxisMax = 45*f0;
ticks = linspace(f0,xAxisMax,5);
labels = cell(length(ticks),1);
labels{1} = 'f_0';
for i=2:(length(ticks)-1)
   labels{i} = sprintf('%df_0',round(ticks(i)/f0));
end
if xAxisMax == Nmax
    labels{length(ticks)} = 'N_{max}';
else
    labels{length(ticks)} = sprintf('%df_0',round(xAxisMax/f0));
end
xticks(ticks)
xticklabels(labels)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Check the variance of u,v,w,zeta as a function of depth
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uvVariance = zeros(length(z),1);
zetaVariance = zeros(length(z),1);
wVariance = zeros(length(z),1);
timeEnsemble = 1:2:(length(t)-1); % 6:min(400,length(t)-1); % Irrelevant if there are no external modes.
for iTime = timeEnsemble
    u = double(squeeze(ncread(file, 'u', [1 1 1 iTime], [length(x) length(y) length(z) 1], [1 1 1 1])));
    v = double(squeeze(ncread(file, 'v', [1 1 1 iTime], [length(x) length(y) length(z) 1], [1 1 1 1])));
    w = double(squeeze(ncread(file, 'w', [1 1 1 iTime], [length(x) length(y) length(z) 1], [1 1 1 1])));
    zeta = double(squeeze(ncread(file, 'zeta', [1 1 1 iTime], [length(x) length(y) length(z) 1], [1 1 1 1])));

    uvVariance = uvVariance + squeeze(mean(mean(u.*u + v.*v,1),2));
    zetaVariance = zetaVariance + squeeze(mean(mean(zeta.*zeta,1),2));
    wVariance = wVariance + squeeze(mean(mean(w.*w,1),2));
end
uvVariance = uvVariance/length(timeEnsemble);
zetaVariance = zetaVariance/length(timeEnsemble);
wVariance = wVariance/length(timeEnsemble);

if shouldWKBScale == 1
    uvVariance = (N0_GM./N) .* uvVariance;
    zetaVariance = (N/N0_GM) .* zetaVariance;
end

% ExpectedUV = E*Phi*( 3*N2/2 - f0*f0 - (B0/2)*f0*sqrt(N2-f0*f0) )/(N2-f0*f0);
% ExpectedZeta = E*Gamma*(1/2 - (B0/2)*(f0/N2)*sqrt(N2-f0*f0))/(N2-f0*f0);
% ExpectedW = E*Gamma*f0*f0*((B0/f0)*sqrt(N2-f0*f0)-1)/(N2-f0*f0);
% 
% ExpectedTotal = ExpectedUV + ExpectedW + N2*ExpectedZeta;

totalVariance = uvVariance + wVariance + N2.*zetaVariance;

figure
subplot(1,4,1)
plot([44 44], [z(1) z(end)], 'k' ,'LineWidth', 2), hold on
plot(1e4*uvVariance,z,'LineWidth', 2)
% plot(1e4*ExpectedUV,z,'LineWidth', 2)
xlabel('cm^2/s^2'), ylabel('depth (m)')
title('HKE')
% title(sprintf('HKE (%.2fGM)',trapz(z,uvVariance)/trapz(z,ExpectedUV) ))

subplot(1,4,2)
plot([44 44], [z(1) z(end)], 'k' ,'LineWidth', 2), hold on
plot(1e4*wVariance,z,'LineWidth', 2)
% plot(1e4*ExpectedW,z,'LineWidth', 2)
xlabel('cm^2/s^2'), set(gca,'YTickLabel',[]);
title('VKE')
% title(sprintf('VKE (%.2fGM)',trapz(z,wVariance)/trapz(z,ExpectedW) ))

subplot(1,4,3)
plot([53 53], [z(1) z(end)], 'k' ,'LineWidth', 2), hold on
plot(zetaVariance,z,'LineWidth', 2), hold on
% plot(ExpectedZeta,z,'LineWidth', 2)
xlabel('m^2'), set(gca,'YTickLabel',[]);
title('isopycnal variance')
% title(sprintf('isopycnal variance (%.2fGM)',trapz(z,zetaVariance)/trapz(z,ExpectedZeta) ))

subplot(1,4,4)
plot([30 30], [z(1) z(end)], 'k' ,'LineWidth', 2), hold on
plot(1e4*totalVariance ,z,'LineWidth', 2), hold on
% plot(1e4*ExpectedTotal ,z,'LineWidth', 2)
xlabel('cm^2 s^{-2}'), set(gca,'YTickLabel',[]);
title('total energy')
% title(sprintf('total energy (%.2fGM)',trapz(z,totalVariance)/trapz(z,ExpectedTotal) ))
