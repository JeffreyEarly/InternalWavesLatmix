file = '/Volumes/home/jearly/InternalWavesLatmixStrained_256_256_80_GM_0.031.nc';
file = '/Volumes/Data/InternalWavesLatmix_128_128_50_GM_0.031_good.nc';

addpath('/Users/jearly/Documents/LatMix/drifters/support')

BoxSize = 1.5e3;
depthCutoff = -50;

xFloat=double(ncread(file, 'x-float'));
yFloat=double(ncread(file, 'y-float'));
zFloat=double(ncread(file, 'z-float'));

t = ncread(file, 'time');
lat0 = ncreadatt(file, '/', 'latitude');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read in the isopycnal following drifters.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xp=double(ncread(file, 'x-position',[1 1 1 1], [Inf Inf Inf length(t)], [1 1 1 1]));
yp=double(ncread(file, 'y-position',[1 1 1 1], [Inf Inf Inf length(t)], [1 1 1 1]));
zp=double(ncread(file, 'z-position',[1 1 1 1], [Inf Inf Inf length(t)], [1 1 1 1]));

zIdx = zeros(size(zp));
zIdx(:,:,2,:)=2;

x = reshape(xp,[],length(t));
y = reshape(yp,[],length(t));
z = reshape(zp,[],length(t));
zIdx = reshape(zIdx,[],length(t));

badDrifters = zeros(size(z));
badDrifters(find(z<depthCutoff))=1;
badDrifters = sum(badDrifters,2);

indices = find( abs(x(:,1)) < BoxSize &  abs(y(:,1)) < BoxSize & zIdx(:,1) == 2 & ~badDrifters);

x = (x(indices,:))';
y = (y(indices,:))';
z = (z(indices,:))';

figure
subplot(2,2,1)
plot(x,y)
% figure
% plot(x,z)
% figure
% scatter(x(1,:),y(1,:))

[minD, maxD, theta] = SecondMomentMatrix( x, y, 'eigen' );
D2 = (minD+maxD)/2;
[p,S,mu]=polyfit(t,D2,1);
kappa_fit = 0.5*p(1)/mu(2);
fprintf('diffusive linear fit: kappa = %f\n', kappa_fit)

subplot(2,2,3)
plot(t, D2)

save('isopycnal_floats_grid.mat','x','y', 'z', 't')

%save(sprintf('%s/isopycnal_floats.mat',FloatFolder),'lat0','t','x','y','z')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read in the isopycnal following, diffusive drifters.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xp=double(ncread(file, 'x-position-diffusive',[1 1 1 1], [Inf Inf Inf length(t)], [1 1 1 1]));
yp=double(ncread(file, 'y-position-diffusive',[1 1 1 1], [Inf Inf Inf length(t)], [1 1 1 1]));
zp=double(ncread(file, 'z-position-diffusive',[1 1 1 1], [Inf Inf Inf length(t)], [1 1 1 1]));

zIdx = zeros(size(zp));
zIdx(:,:,2,:)=2;

x = reshape(xp,[],length(t));
y = reshape(yp,[],length(t));
z = reshape(zp,[],length(t));
zIdx = reshape(zIdx,[],length(t));

badDrifters = zeros(size(z));
badDrifters(find(z<depthCutoff | z > 0))=1;
badDrifters = sum(badDrifters,2);

indices = find( abs(x(:,1)) < BoxSize &  abs(y(:,1)) < BoxSize & zIdx(:,1) == 2 & ~badDrifters);

x = (x(indices,:))';
y = (y(indices,:))';
z = (z(indices,:))';


subplot(2,2,2)
plot(x,y)
% figure
% plot(x,z)
% figure
% scatter(x(1,:),y(1,:))

[minD, maxD, theta] = SecondMomentMatrix( x, y, 'eigen' );
D2 = (minD+maxD)/2;
[p,S,mu]=polyfit(t,D2,1);
kappa_fit = 0.5*p(1)/mu(2);
fprintf('diffusive linear fit: kappa = %f\n', kappa_fit)
subplot(2,2,4), plot(t, D2)

save('isopycnal_diffusive_floats_grid.mat','x','y', 'z', 't')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read in the fixed depth floats.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stride = 4;
xp=double(ncread(file, 'x-position-fixed-depth',[1 1 1 1], [length(xFloat)/stride length(yFloat)/stride Inf length(t)], [stride stride 1 1]));
yp=double(ncread(file, 'y-position-fixed-depth',[1 1 1 1], [length(xFloat)/stride length(yFloat)/stride Inf length(t)], [stride stride 1 1]));
zp=double(ncread(file, 'z-position-fixed-depth',[1 1 1 1], [length(xFloat)/stride length(yFloat)/stride Inf length(t)], [stride stride 1 1]));

zIdx = zeros(size(zp));
zIdx(:,:,2,:)=2;

x = reshape(xp,[],length(t));
y = reshape(yp,[],length(t));
z = reshape(zp,[],length(t));
zIdx = reshape(zIdx,[],length(t));

badDrifters = zeros(size(z));
badDrifters(find(z<-50))=1;
badDrifters = sum(badDrifters,2);

indices = find( abs(x(:,1)) < BoxSize &  abs(y(:,1)) < BoxSize & zIdx(:,1) == 2 & ~badDrifters);

x = (x(indices,:))';
y = (y(indices,:))';
z = (z(indices,:))';

figure
subplot(2,2,1)
plot(x,y)
% figure
% plot(x,z)
% figure
% scatter(x(1,:),z(1,:))

[minD, maxD, theta] = SecondMomentMatrix( x, y, 'eigen' );
D2 = (minD+maxD)/2;
[p,S,mu]=polyfit(t,D2,1);
kappa_fit = 0.5*p(1)/mu(2);
fprintf('diffusive linear fit: kappa = %f\n', kappa_fit)

subplot(2,2,3), plot(t, D2)

save('fixed_depth_floats_grid.mat','x','y', 'z', 't')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read in the drogued drifters.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xp=double(ncread(file, 'x-position-drifter',[1 1 1], [length(xFloat)/stride length(yFloat)/stride length(t)], [stride stride 1]));
yp=double(ncread(file, 'y-position-drifter',[1 1 1], [length(xFloat)/stride length(yFloat)/stride length(t)], [stride stride 1]));

x = reshape(xp,[],length(t));
y = reshape(yp,[],length(t));

indices = find( abs(x(:,1)) < BoxSize &  abs(y(:,1)) < BoxSize );

x = (x(indices,:))';
y = (y(indices,:))';

subplot(2,2,2)
plot(x,y)

[minD, maxD, theta] = SecondMomentMatrix( x, y, 'eigen' );
D2 = (minD+maxD)/2;
[p,S,mu]=polyfit(t,D2,1);
kappa_fit = 0.5*p(1)/mu(2);
fprintf('diffusive linear fit: kappa = %f\n', kappa_fit)
subplot(2,2,4), plot(t, D2)

save('drifters_grid.mat','x','y','t')