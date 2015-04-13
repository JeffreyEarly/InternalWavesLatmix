file = '/Volumes/home/jearly/InternalWavesLatmixStrained_256_256_80_GM_0.031.nc';
FloatFolder = '/Users/jearly/Documents/LatMix/model/Square_32GM_Strained_HighResolution';

if exist(FloatFolder,'dir') == 0
	mkdir(FloatFolder);
end

xFloat=double(ncread(file, 'x-float'));
yFloat=double(ncread(file, 'y-float'));
zFloat=double(ncread(file, 'z-float'));

t = ncread(file, 'time');
lat0 = ncreadatt(file, '/', 'latitude');

xIdx = [31 31 29 30 31 32 33 31 31];
yIdx = [29 30 31 31 31 31 31 32 33];
zIdx = [2 2 2 2 2 2 2 2 2];

xIdx = 2*[31 31 29 30 31 32 33 31 31]-30;
yIdx = 2*[29 30 31 31 31 31 31 32 33]-30;
zIdx = [2 2 2 2 2 2 2 2 2];

xIdx = 4*[3 3 1 2 3 4 5 3 3]; 
yIdx = 4*[1 2 3 3 3 3 3 4 5]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read in the isopycnal following drifters.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xp=double(ncread(file, 'x-position',[1 1 1 1], [Inf Inf Inf length(t)], [1 1 1 1]));
yp=double(ncread(file, 'y-position',[1 1 1 1], [Inf Inf Inf length(t)], [1 1 1 1]));
zp=double(ncread(file, 'z-position',[1 1 1 1], [Inf Inf Inf length(t)], [1 1 1 1]));

x = zeros(length(t),length(xIdx));
y = zeros(length(t),length(xIdx));
z = zeros(length(t),length(xIdx));

for i=1:length(xIdx)
   x(:,i) = xp(xIdx(i),yIdx(i),zIdx(i),:);
   y(:,i) = yp(xIdx(i),yIdx(i),zIdx(i),:);
   z(:,i) = zp(xIdx(i),yIdx(i),zIdx(i),:);
end

figure
subplot(2,2,1)
plot(x,y)
save(sprintf('%s/isopycnal_floats.mat',FloatFolder),'lat0','t','x','y','z')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read in the isopycnal following, diffusive drifters.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xp=double(ncread(file, 'x-position-diffusive',[1 1 1 1], [Inf Inf Inf length(t)], [1 1 1 1]));
yp=double(ncread(file, 'y-position-diffusive',[1 1 1 1], [Inf Inf Inf length(t)], [1 1 1 1]));
zp=double(ncread(file, 'z-position-diffusive',[1 1 1 1], [Inf Inf Inf length(t)], [1 1 1 1]));

x = zeros(length(t),length(xIdx));
y = zeros(length(t),length(xIdx));
z = zeros(length(t),length(xIdx));

for i=1:length(xIdx)
   x(:,i) = xp(xIdx(i),yIdx(i),zIdx(i),:);
   y(:,i) = yp(xIdx(i),yIdx(i),zIdx(i),:);
   z(:,i) = zp(xIdx(i),yIdx(i),zIdx(i),:);
end

subplot(2,2,2)
plot(x,y)
save(sprintf('%s/isopycnal_diffusive_floats.mat',FloatFolder),'lat0','t','x','y','z')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read in the fixed depth floats.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xp=double(ncread(file, 'x-position-fixed-depth',[1 1 1 1], [Inf Inf Inf length(t)], [1 1 1 1]));
yp=double(ncread(file, 'y-position-fixed-depth',[1 1 1 1], [Inf Inf Inf length(t)], [1 1 1 1]));
zp=double(ncread(file, 'z-position-fixed-depth',[1 1 1 1], [Inf Inf Inf length(t)], [1 1 1 1]));

x = zeros(length(t),length(xIdx));
y = zeros(length(t),length(xIdx));
z = zeros(length(t),length(xIdx));

for i=1:length(xIdx)
   x(:,i) = xp(xIdx(i),yIdx(i),zIdx(i),:);
   y(:,i) = yp(xIdx(i),yIdx(i),zIdx(i),:);
   z(:,i) = zp(xIdx(i),yIdx(i),zIdx(i),:);
end

subplot(2,2,3)
plot(x,y)
save(sprintf('%s/fixed_depth_floats.mat',FloatFolder),'lat0','t','x','y','z')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read in the drogued drifters.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xp=double(ncread(file, 'x-position-drifter',[1 1 1], [Inf Inf length(t)], [1 1 1]));
yp=double(ncread(file, 'y-position-drifter',[1 1 1], [Inf Inf length(t)], [1 1 1]));

x = zeros(length(t),length(xIdx));
y = zeros(length(t),length(xIdx));

for i=1:length(xIdx)
   x(:,i) = xp(xIdx(i),yIdx(i),:);
   y(:,i) = yp(xIdx(i),yIdx(i),:);
end

subplot(2,2,4)
plot(x,y)
save(sprintf('%s/drifters.mat',FloatFolder),'lat0','t','x','y')