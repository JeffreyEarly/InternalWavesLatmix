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

N = 64; % probably should be 128
aspectRatio = 4;

L = 25e3;
Lx = aspectRatio*L;
Ly = L;
Lz = 5000;

Nx = aspectRatio*N;
Ny = N;
Nz = 32; 
nModes = 64;

latitude = 31;
GMReferenceLevel = 1.0;

kappa = 5e-6;
outputInterval = 15*60;
% maxTime = 6.0*86400;
maxTime = 10*outputInterval;
interpolationMethod = 'spline';

shouldOutputEulerianFields = 1;
shouldOutputFloats = 1;
shouldOutputDrifters = 1;

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

[rho, N2, zIn] = InternalModes.StratificationProfileWithName('latmix-site1');

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
    wavemodel = InternalWaveModelArbitraryStratification([Lx, Ly, Lz], [Nx, Ny], rho, z, latitude, nModes,  'method', 'wkbAdaptiveSpectral');
end

wavemodel.FillOutWaveSpectrum();
wavemodel.InitializeWithGMSpectrum(GMReferenceLevel,1);
wavemodel.ShowDiagnostics();

period = 2*pi/wavemodel.Nmax;
if shouldOutputEulerianFields == 1
    [u,v,w] = wavemodel.VelocityFieldAtTime(0.0);
    U = max(max(max( sqrt(u.*u + v.*v) )));
    W = max(max(max( abs(w) )));
else
    U = 0.1;
end
fprintf('Max fluid velocity: U=%.2f cm/s, W=%.2f\n',U*100,W*100);


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

zIsopycnal = wavemodel.PlaceParticlesOnIsopycnal(x_float,y_float,z_float,interpolationMethod,1e-6);

y0 = cat(2, cat(1,x_float,x_float), cat(1,y_float,y_float),cat(1,z_float,zIsopycnal));
drifterIndices = (1:length(z_float))';
floatIndices = drifterIndices+length(z_float);
nDrifters = numel(drifterIndices);
nFloats = numel(floatIndices);
shouldUseW = cat(1,0*z_float,zIsopycnal);

f = @(t,y) wavemodel.VelocityAtTimePositionVector(t,y,interpolationMethod,shouldUseW);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Determine the proper time interval
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfl = 0.25;
advectiveDT = cfl*(wavemodel.x(2)-wavemodel.x(1))/U;
advectiveWDT = cfl*(min(abs(diff(z))))/W;
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

filepath = sprintf('%s/LatmixSite1Experiment_%s_%dx%dx%d.nc', outputfolder,datestr(datetime('now'),'yyyy-mm-ddTHHMMSS'),Nx,Ny,Nz);

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
if shouldOutputFloats == 1
    floatDimID = netcdf.defDim(ncid, 'float_id', nFloats);
    xFloatID = netcdf.defVar(ncid, 'x-position-float', ncPrecision, [floatDimID,tDimID]);
    yFloatID = netcdf.defVar(ncid, 'y-position-float', ncPrecision, [floatDimID,tDimID]);
    zFloatID = netcdf.defVar(ncid, 'z-position-float', ncPrecision, [floatDimID,tDimID]);
    netcdf.putAtt(ncid,xFloatID, 'units', 'm');
    netcdf.putAtt(ncid,yFloatID, 'units', 'm');
    netcdf.putAtt(ncid,zFloatID, 'units', 'm');
end

% Define the *drifter* dimensions
if shouldOutputDrifters == 1
    drifterDimID = netcdf.defDim(ncid, 'drifter_id', nDrifters);
    xDrifterID = netcdf.defVar(ncid, 'x-position-drifter', ncPrecision, [drifterDimID,tDimID]);
    yDrifterID = netcdf.defVar(ncid, 'y-position-drifter', ncPrecision, [drifterDimID,tDimID]);
    zDrifterID = netcdf.defVar(ncid, 'z-position-drifter', ncPrecision, [drifterDimID,tDimID]);
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
integrator = Integrator( f, y0, deltaT);
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
        [u,v,w,zeta,rho]=wavemodel.VariableFieldsAtTime(t(iTime),'u', 'v', 'w', 'zeta', 'rho_prime');
        
        netcdf.putVar(ncid, uVarID, [0 0 0 iTime-1], [wavemodel.Nx wavemodel.Ny wavemodel.Nz 1], u);
        netcdf.putVar(ncid, vVarID, [0 0 0 iTime-1], [wavemodel.Nx wavemodel.Ny wavemodel.Nz 1], v);
        netcdf.putVar(ncid, wVarID, [0 0 0 iTime-1], [wavemodel.Nx wavemodel.Ny wavemodel.Nz 1], w);
        netcdf.putVar(ncid, zetaVarID, [0 0 0 iTime-1], [wavemodel.Nx wavemodel.Ny wavemodel.Nz 1], zeta);
        netcdf.putVar(ncid, rhoVarID, [0 0 0 iTime-1], [wavemodel.Nx wavemodel.Ny wavemodel.Nz 1], rho);
    end
    netcdf.putVar(ncid, tVarID, iTime-1, 1, t(iTime));
    
    if shouldOutputDrifters == 1 && shouldOutputFloats == 1
        p = integrator.StepForwardToTime(t(iTime));
        
        x_drifter = p(drifterIndices,1);
        y_drifter = p(drifterIndices,2);
        z_drifter = p(drifterIndices,3);
        
        x_float = p(floatIndices,1);
        y_float = p(floatIndices,2);
        z_float = p(floatIndices,3);
        
        netcdf.putVar(ncid, setprecision(xDrifterID), [0 iTime-1], [nDrifters 1], x_drifter);
        netcdf.putVar(ncid, setprecision(yDrifterID), [0 iTime-1], [nDrifters 1], y_drifter);
        netcdf.putVar(ncid, setprecision(zDrifterID), [0 iTime-1], [nDrifters 1], z_drifter);
        
        netcdf.putVar(ncid, setprecision(xFloatID), [0 iTime-1], [nFloats 1], x_float);
        netcdf.putVar(ncid, setprecision(yFloatID), [0 iTime-1], [nFloats 1], y_float);
        netcdf.putVar(ncid, setprecision(zFloatID), [0 iTime-1], [nFloats 1], z_float);
    end

end
% profile viewer
fprintf('Ending numerical simulation on %s\n', datestr(datetime('now')));

netcdf.close(ncid);	
