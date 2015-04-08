file = '/Users/jearly/Desktop/InternalWavesLatmix_128_128_50_GM_0.013.nc';
file = '/Volumes/Data/InternalWavesLatmix_256_256_50_GM_0.062.nc';
file = '/Volumes/jearly/Desktop/InternalWavesLatmix_256_256_50_GM_0.062.nc';
FloatFolder = '/Users/jearly/Desktop/InternalWavesLatmix_128_128_50_GM_0.013';
FloatFolder = '/Users/jearly/Documents/LatMix/model/eighthGM';
FloatFolder = '/Users/jearly/Documents/LatMix/model/sixteenthGM';
%FloatFolder = '/Users/jearly/Documents/LatMix/model/thirtysecondGM';

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

xIdx = 2*[31 31 29 30 31 32 33 31 31];
yIdx = 2*[29 30 31 31 31 31 31 32 33];
zIdx = [2 2 2 2 2 2 2 2 2];

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

x = zeros(length(t),length(xIdx));
y = zeros(length(t),length(xIdx));
z = zeros(length(t),length(xIdx));

for i=1:length(xIdx)
   x(:,i) = xp(xIdx(i),yIdx(i),zIdx(i),:);
   y(:,i) = yp(xIdx(i),yIdx(i),zIdx(i),:);
   z(:,i) = zp(xIdx(i),yIdx(i),zIdx(i),:);
end

xp=double(ncread(file, 'x-position-diffusive',[1 1 1 1], [Inf Inf Inf length(t)], [1 1 1 1]));
yp=double(ncread(file, 'y-position-diffusive',[1 1 1 1], [Inf Inf Inf length(t)], [1 1 1 1]));
zp=double(ncread(file, 'z-position-diffusive',[1 1 1 1], [Inf Inf Inf length(t)], [1 1 1 1]));

subplot(2,2,2)
plot(x,y)
save(sprintf('%s/isopycnal_diffusive_floats.mat',FloatFolder),'lat0','t','x','y','z')

x = zeros(length(t),length(xIdx));
y = zeros(length(t),length(xIdx));
z = zeros(length(t),length(xIdx));

for i=1:length(xIdx)
   x(:,i) = xp(xIdx(i),yIdx(i),zIdx(i),:);
   y(:,i) = yp(xIdx(i),yIdx(i),zIdx(i),:);
   z(:,i) = zp(xIdx(i),yIdx(i),zIdx(i),:);
end

xp=double(ncread(file, 'x-position-fixed-depth',[1 1 1 1], [Inf Inf Inf length(t)], [1 1 1 1]));
yp=double(ncread(file, 'y-position-fixed-depth',[1 1 1 1], [Inf Inf Inf length(t)], [1 1 1 1]));
zp=double(ncread(file, 'z-position-fixed-depth',[1 1 1 1], [Inf Inf Inf length(t)], [1 1 1 1]));

subplot(2,2,3)
plot(x,y)
save(sprintf('%s/fixed_depth_floats.mat',FloatFolder),'lat0','t','x','y','z')

x = zeros(length(t),length(xIdx));
y = zeros(length(t),length(xIdx));
z = zeros(length(t),length(xIdx));

for i=1:length(xIdx)
   x(:,i) = xp(xIdx(i),yIdx(i),zIdx(i),:);
   y(:,i) = yp(xIdx(i),yIdx(i),zIdx(i),:);
   z(:,i) = zp(xIdx(i),yIdx(i),zIdx(i),:);
end

xp=double(ncread(file, 'x-position-drifter',[1 1 1 1], [Inf Inf Inf length(t)], [1 1 1 1]));
yp=double(ncread(file, 'y-position-drifter',[1 1 1 1], [Inf Inf Inf length(t)], [1 1 1 1]));
zp=double(ncread(file, 'z-position-drifter',[1 1 1 1], [Inf Inf Inf length(t)], [1 1 1 1]));

subplot(2,2,4)
plot(x,y)
save(sprintf('%s/drifters.mat',FloatFolder),'lat0','t','x','y','z')