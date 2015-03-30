file = '/Users/jearly/Desktop/InternalWavesLatmix_128_128_50_GM_0.013.nc';
FloatFolder = '/Users/jearly/Desktop/InternalWavesLatmix_128_128_50_GM_0.013';
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

xp=double(ncread(file, 'x-position'));
yp=double(ncread(file, 'y-position'));
zp=double(ncread(file, 'z-position'));

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

xp=double(ncread(file, 'x-position-diffusive'));
yp=double(ncread(file, 'y-position-diffusive'));
zp=double(ncread(file, 'z-position-diffusive'));

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

xp=double(ncread(file, 'x-position-fixed-depth'));
yp=double(ncread(file, 'y-position-fixed-depth'));
zp=double(ncread(file, 'z-position-fixed-depth'));

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

xp=double(ncread(file, 'x-position-drifter'));
yp=double(ncread(file, 'y-position-drifter'));
zp=double(ncread(file, 'z-position-drifter'));

subplot(2,2,4)
plot(x,y)
save(sprintf('%s/drifters.mat',FloatFolder),'lat0','t','x','y','z')