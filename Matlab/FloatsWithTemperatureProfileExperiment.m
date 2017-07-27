%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DiffusivityExperiment
%
% This script uses the InternalWaveModel to generate the time-evolution of
% a Garrett-Munk spectrum of internal waves and save the output to a NetCDF
% file. Advects particles as well.
%
% Jeffrey J. Early
% jeffrey@jeffreyearly.com
%
% April 12th, 2017      Version 1.0

N = 128; % probably should be 128
aspectRatio = 4;

L = 25e3;
Lx = aspectRatio*L;
Ly = L;
Lz = 5000;

Nx = aspectRatio*N;
Ny = N;
Nz = 64; 
nModes = 64;
nEVP = 512; % probably should be 512
nGrid = 2^14+1;

latitude = 31;
GMReferenceLevel = 1.3;

kappa = 5e-6;
outputInterval = 15*60;
maxTime = 6.0*86400; %10*outputInterval;
interpolationMethod = 'spline';

shouldOutputEulerianFields = 1;
shouldOutputFloats = 0;
shouldOutputDiffusiveFloats = 0;
shouldOutputDrifters = 0;

outputfolder = '/Volumes/OceanTransfer';
% outputfolder = '/Users/jearly/Desktop';

precision = 'single';

if strcmp(precision,'single')
    ncPrecision = 'NC_FLOAT';
    setprecision = @(x) single(x);
    bytePerFloat = 4;
else
    ncPrecision = 'NC_DOUBLE';
    setprecision = @(x) double(x);
    bytePerFloat = 8;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Build the Latmix stratification profile
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Notes for this profile are here:
% /Users/jearly/Dropbox/Documents/Notes/latmix-glider-profiles/ADifferentAnalyticalProfile/ADifferentAnalyticalProfile.tex

rho0 = 1025;
g = 9.81;

delta_p = 0.9;
L_s = 8;
L_d = 100;
z_p = -17;
D = Lz; % intended to be 5000m

N0 = 2.8e-3; % Surface
Nq = 1.4e-2; % Pycnocline portion to exponentials
Np = 4.7e-2;
Nd = 1.75e-3;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Now let's create a stretched coordinate and see how well it does
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxDepth = -50;
z_lin = linspace(0,maxDepth,1e5);

N_scaled = sqrt(N2(z_lin)/N2(0));
s = cumtrapz( z_lin, N_scaled );
sGrid = linspace(min(s), max(s), Nz);
z = interp1(s, z_lin, sGrid );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shouldUseGMSpectrum = 1;

if ~exist('wavemodel','var')
    wavemodel = InternalWaveModelArbitraryStratification([Lx, Ly, Lz], [Nx, Ny], rho, z, latitude, nModes,  'method', 'wkbSpectral', 'nEVP', nEVP, 'nGrid', nGrid);
end

wavemodel.FillOutWaveSpectrum();
wavemodel.InitializeWithGMSpectrum(GMReferenceLevel,0);
wavemodel.ShowDiagnostics();

period = 2*pi/wavemodel.Nmax;
if shouldOutputEulerianFields == 1
    [u,v] = wavemodel.VelocityFieldAtTime(0.0);
    U = max(max(max( sqrt(u.*u + v.*v) )));
else
    U = 0.1;
end
fprintf('Max fluid velocity: %.2f cm/s\n',U*100);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create floats/drifters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dx = wavemodel.x(2)-wavemodel.x(1);
dy = wavemodel.y(2)-wavemodel.y(1);
nLevels = 5;
N = floor(N/3);
x_float = (0:N-1)*dx;
y_float = (0:N-1)*dy;
z_float = -30;

% nudge towards the center of the domain. This isn't necessary, but does
% prevent the spline interpolation from having to worry about the
% boundaries.
x_float = x_float + (max(wavemodel.x) - max(x_float))/2;
y_float = y_float + (max(wavemodel.y) - max(y_float))/2;

[x_float,y_float,z_float] = ndgrid(x_float,y_float,z_float);
x_float = reshape(x_float,[],1);
y_float = reshape(y_float,[],1);
z_float = reshape(z_float,[],1);
nFloats = numel(x_float);

ymin = [-Inf -Inf];
ymax = [Inf Inf];
kappa_vector = [0 0];
p0 = cat(2, x_float, y_float);

f = @(t,y) FluxForDrifter(t,y,z_float,wavemodel, interpolationMethod);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Determine the proper time interval
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfl = 0.25;
advectiveDT = cfl*(wavemodel.x(2)-wavemodel.x(1))/U;
oscillatoryDT = period/8;
if advectiveDT < oscillatoryDT
    fprintf('Using the advective dt: %.2f\n',advectiveDT);
    deltaT = advectiveDT;
else
    fprintf('Using the oscillatory dt: %.2f\n',oscillatoryDT);
    deltaT = oscillatoryDT;
end

deltaT = outputInterval/ceil(outputInterval/deltaT);
fprintf('Rounding to match the output interval dt: %.2f\n',deltaT);

t = (0:outputInterval:maxTime)';
if t(end) < period
    t(end+1) = period;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create a NetCDF file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filepath = sprintf('%s/FloatsWithTemperatureProfileExperiment_%s_%dx%dx%d.nc', outputfolder,datestr(datetime('now'),'yyyy-mm-ddTHHMMSS'),Nx,Ny,Nz);

% Apple uses 1e9 bytes as 1 GB (not the usual multiples of 2 definition)
totalFields = 4;
totalSize = totalFields*bytePerFloat*length(t)*(wavemodel.Nx)*(wavemodel.Ny)*(wavemodel.Nz)/1e9;
fprintf('Writing output file to %s\nExpected file size is %.2f GB.\n',filepath,totalSize);

cmode = netcdf.getConstant('CLOBBER');
cmode = bitor(cmode,netcdf.getConstant('SHARE'));
ncid = netcdf.create(filepath, cmode);

% Define the dimensions
xDimID = netcdf.defDim(ncid, 'x', wavemodel.Nx);
yDimID = netcdf.defDim(ncid, 'y', wavemodel.Ny);
zDimID = netcdf.defDim(ncid, 'z', wavemodel.Nz);
tDimID = netcdf.defDim(ncid, 't', netcdf.getConstant('NC_UNLIMITED'));

% Define the coordinate variables
xVarID = netcdf.defVar(ncid, 'x', ncPrecision, xDimID);
yVarID = netcdf.defVar(ncid, 'y', ncPrecision, yDimID);
zVarID = netcdf.defVar(ncid, 'z', ncPrecision, zDimID);
tVarID = netcdf.defVar(ncid, 't', ncPrecision, tDimID);
netcdf.putAtt(ncid,xVarID, 'units', 'm');
netcdf.putAtt(ncid,yVarID, 'units', 'm');
netcdf.putAtt(ncid,zVarID, 'units', 'm');
netcdf.putAtt(ncid,tVarID, 'units', 's');

% Define the dynamical variables
if shouldOutputEulerianFields == 1
    uVarID = netcdf.defVar(ncid, 'u', ncPrecision, [xDimID,yDimID,zDimID,tDimID]);
    vVarID = netcdf.defVar(ncid, 'v', ncPrecision, [xDimID,yDimID,zDimID,tDimID]);
    wVarID = netcdf.defVar(ncid, 'w', ncPrecision, [xDimID,yDimID,zDimID,tDimID]);
    zetaVarID = netcdf.defVar(ncid, 'zeta', ncPrecision, [xDimID,yDimID,zDimID,tDimID]);
    rhoVarID = netcdf.defVar(ncid, 'rho', ncPrecision, [xDimID,yDimID,zDimID,tDimID]);
    netcdf.putAtt(ncid,uVarID, 'units', 'm/s');
    netcdf.putAtt(ncid,vVarID, 'units', 'm/s');
    netcdf.putAtt(ncid,wVarID, 'units', 'm/s');
    netcdf.putAtt(ncid,zetaVarID, 'units', 'm');
    netcdf.putAtt(ncid,rhoVarID, 'units', 'kg/m^3');
end

% Define the *float* dimensions
if shouldOutputDrifters == 1
    xDrifterID = netcdf.defVar(ncid, 'x-position-drifter', ncPrecision, [floatDimID,tDimID]);
    yDrifterID = netcdf.defVar(ncid, 'y-position-drifter', ncPrecision, [floatDimID,tDimID]);
    zDrifterID = netcdf.defVar(ncid, 'z-position-drifter', ncPrecision, [floatDimID,tDimID]);
    densityDrifterID = netcdf.defVar(ncid, 'density-drifter', ncPrecision, [floatDimID,tDimID]);
    netcdf.putAtt(ncid,xDrifterID, 'units', 'm');
    netcdf.putAtt(ncid,yDrifterID, 'units', 'm');
    netcdf.putAtt(ncid,zDrifterID, 'units', 'm');
end

% Write the density profile
n2VarID = netcdf.defVar(ncid, 'N2', ncPrecision, zDimID);
netcdf.putAtt(ncid,n2VarID, 'units', '1/s^2');

% Write some metadata
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'latitude', latitude);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'GMReferenceLevel', GMReferenceLevel);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'Model', 'Created from InternalWaveModel.m written by Jeffrey J. Early.');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'ModelVersion', wavemodel.version);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'CreationDate', datestr(datetime('now')));

netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'kappa', kappa);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'interpolation-method', interpolationMethod);

% End definition mode
netcdf.endDef(ncid);

% Add the data for the coordinate variables
netcdf.putVar(ncid, setprecision(xVarID), wavemodel.x);
netcdf.putVar(ncid, setprecision(yVarID), wavemodel.y);
netcdf.putVar(ncid, setprecision(zVarID), wavemodel.z);
netcdf.putVar(ncid, n2VarID, wavemodel.N2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Run the model, and write the output to NetCDF
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
startTime = datetime('now');
fprintf('Starting numerical simulation on %s\n', datestr(startTime));
integrator = IntegratorWithDiffusivity( f, p0, deltaT, kappa_vector, ymin, ymax);
% profile on
for iTime=1:length(t)
    if iTime == 2
       startTime = datetime('now'); 
    end
    if iTime == 3 || mod(iTime,10) == 0
        timePerStep = (datetime('now')-startTime)/(iTime-2);
        timeRemaining = (length(t)-iTime+1)*timePerStep;   
        fprintf('\twriting values time step %d of %d to file. Estimated finish time %s (%s from now)\n', iTime, length(t), datestr(datetime('now')+timeRemaining), datestr(timeRemaining, 'HH:MM:SS')) ;
    end

    if shouldOutputEulerianFields == 1
        [u,v]=wavemodel.VelocityFieldAtTime(t(iTime));
        [w,zeta] = wavemodel.VerticalFieldsAtTime(t(iTime));
        rho = wavemodel.DensityAtTime(t(iTime));
        
        netcdf.putVar(ncid, uVarID, [0 0 0 iTime-1], [wavemodel.Nx wavemodel.Ny wavemodel.Nz 1], u);
        netcdf.putVar(ncid, vVarID, [0 0 0 iTime-1], [wavemodel.Nx wavemodel.Ny wavemodel.Nz 1], v);
        netcdf.putVar(ncid, wVarID, [0 0 0 iTime-1], [wavemodel.Nx wavemodel.Ny wavemodel.Nz 1], w);
        netcdf.putVar(ncid, zetaVarID, [0 0 0 iTime-1], [wavemodel.Nx wavemodel.Ny wavemodel.Nz 1], zeta);
        netcdf.putVar(ncid, rhoVarID, [0 0 0 iTime-1], [wavemodel.Nx wavemodel.Ny wavemodel.Nz 1], rho);
    end
    netcdf.putVar(ncid, tVarID, iTime-1, 1, t(iTime));
    
    if shouldOutputDrifters == 1
        p = integrator.StepForwardToTime(t(iTime));
        netcdf.putVar(ncid, setprecision(xDrifterID), [0 iTime-1], [nFloats 1], p(:,7));
        netcdf.putVar(ncid, setprecision(yDrifterID), [0 iTime-1], [nFloats 1], p(:,8));
        netcdf.putVar(ncid, setprecision(zDrifterID), [0 iTime-1], [nFloats 1], z_float);
        netcdf.putVar(ncid, setprecision(densityDrifterID), [0 iTime-1], [nFloats 1], wavemodel.DensityAtTimePosition(t(iTime),p(:,7),p(:,8),z_float)-wavemodel.rho0);
    end

end
% profile viewer
fprintf('Ending numerical simulation on %s\n', datestr(datetime('now')));

netcdf.close(ncid);	
