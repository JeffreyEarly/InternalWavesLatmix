//
//  main.m
//  InternalWavesLatmix
//
//  Created by Jeffrey Early on 3/28/15.
//  Copyright (c) 2015 Early Innovations, LLC. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>
#import <GLOceanKit/GLOceanKit.h>

int main(int argc, const char * argv[])
{
    @autoreleasepool {
        GLFloat latitude = 31;
        GLFloat width = 15e3;
        GLFloat height = 15e3;
        NSUInteger Nx = 256;
        NSUInteger Ny = 256;
        NSUInteger Nz_in = 512; // Number of grid points upon which to project the input profile (512 rec.)
        
        NSUInteger Nz_out = 80; // Number of grid points and range for the output
        GLFloat minDepth = -100;
        GLFloat maxDepth = 0;
        
        GLFloat maxWavePeriods = 1.0;
        GLFloat horizontalFloatSpacingInMeters = 125;
        GLFloat sampleTimeInMinutes = 15;
        GLFloat energyLevel = 1./16.;
        
        NSString *initialConditionsFile = [[NSSearchPathForDirectoriesInDomains(NSDesktopDirectory, NSUserDomainMask, YES) objectAtIndex:0] stringByAppendingPathComponent: [NSString stringWithFormat: @"InternalWavesLatmix_%lu_%lu_%lu_%luKM.internalwaves", Nx, Ny, Nz_out,(NSUInteger) (width/1e3)]];
        GLFloat f0 = 2*(7.2921e-5)*sin(latitude*M_PI/180);
        
        /************************************************************************************************/
        /*		Create a density profile and compute the internal wave phases                           */
        /************************************************************************************************/
        GLInternalWaveInitialization *wave;
        GLDimension *xDim, *yDim, *zDim;
        GLEquation *equation;
        
        NSFileManager *manager = [[NSFileManager alloc] init];
        if ([manager fileExistsAtPath: initialConditionsFile isDirectory: nil])
        {
            wave = [NSKeyedUnarchiver unarchiveObjectWithFile: initialConditionsFile];
            if (!wave) {
                NSLog(@"Failed to load wave file.");
                return 0;
            }
            zDim = wave.fullDimensions[0];
            xDim = wave.fullDimensions[1];
            yDim = wave.fullDimensions[2];
            equation = wave.equation;
        }
        else
        {
            equation = [[GLEquation alloc] init];
            
            GLNetCDFFile *profile = [[GLNetCDFFile alloc] initWithURL:[NSURL URLWithString: @"/Users/jearly/Documents/Models/InternalWaves/Latmix2011Site1Profile_AllPoints.nc"] forEquation:equation];
            GLFunction *rho_full = profile.variables[0];
            GLFunction *N2_full = profile.variables[1];
            GLDimension *zDim_full = rho_full.dimensions[0];
            
            // Interpolate the density onto the reduced grid (our z input grid).
            // Do I want to interpolate to a Chebyshev grid?
            GLDimension *zDimRho = [[GLDimension alloc] initDimensionWithGrid: kGLChebyshevEndpointGrid nPoints: Nz_in domainMin: zDim_full.domainMin length: zDim_full.domainLength];
            GLFunction *z = [GLFunction functionOfRealTypeFromDimension:zDimRho withDimensions:@[zDimRho] forEquation:equation];
            GLFunction *rho_bar = [rho_full interpolateAtPoints:@[z]];
            rho_bar.name = @"rho_bar";
            [rho_bar solve];
            
            // Now create a stretched z-output grid
            GLFunction *N_scaled = [[N2_full dividedBy: [N2_full min]] sqrt];
            GLFunction *s = [N_scaled integrate]; // we now have s(z)
            
            // Now determine what s is at the min and max depths
            GLFloat sAtMinDepth = 0.0;
            GLFloat sAtMaxDepth = 0.0;
            [s solve]; // Must solve before trying to use its data to initialize a dimension.
            GLDimension *sDim = [[GLDimension alloc] initWithNPoints: s.nDataPoints values: s.data];
            GLFunction *zOfs = [GLFunction functionOfRealTypeWithDimensions: @[sDim] forEquation: equation];
            for (NSUInteger i=0; i<zOfs.nDataPoints; i++) {
                zOfs.pointerValue[i] = [zDim_full valueAtIndex: i];
                if ([zDim_full valueAtIndex: i] <= minDepth) {
                    sAtMinDepth = s.pointerValue[i];
                }
                if ([zDim_full valueAtIndex: i] >= maxDepth && !sAtMaxDepth) {
                    sAtMaxDepth = s.pointerValue[i];
                }
            }
            
            GLDimension *sDimOut = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: Nz_out domainMin: sAtMinDepth length: sAtMaxDepth-sAtMinDepth];
            GLFunction *sOut = [GLFunction functionOfRealTypeFromDimension: sDimOut withDimensions:@[sDimOut] forEquation:equation];
            GLFunction *zInterp = [zOfs interpolateAtPoints:@[sOut]];
            
            [zInterp solve]; // required!!!
            zDim = [[GLDimension alloc] initWithNPoints: zInterp.nDataPoints values: zInterp.data];
            zDim.name = @"z";
            xDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: Nx domainMin: -width/2 length: width];
            xDim.name = @"x";
            yDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: Ny domainMin: -height/2 length: height];
            yDim.name = @"y";
            
            // The ordering the dimensions is very deliberate here, for two reasons:
            // 1. The z-dimension is in the first position, so that the horizontal fft will act on contiguous chunks of memory, and
            // 2. The last two dimensions are ordered (x,y) to appease pcolor, meshgrid, and all the standard matlab formating.
            wave = [[GLInternalWaveInitialization alloc] initWithDensityProfile: rho_bar fullDimensions:@[zDim, xDim, yDim] latitude:latitude equation:equation];
            
            if (![NSKeyedArchiver archiveRootObject: wave toFile: initialConditionsFile]) {
                NSLog(@"Failed to save restart file.");
            }
        }
        
        [wave createGarrettMunkSpectrumWithEnergy: energyLevel];
        
        // Create the time dimension
        GLFloat maxTime = maxWavePeriods*2*M_PI/f0;
        GLFloat sampleTime = sampleTimeInMinutes*60;
        NSLog(@"Maximum wave period %d @ %02d:%02d (HH:MM)", (int) floor(maxTime/86400), ((int) floor(maxTime/3600))%24, ((int) floor(maxTime/60))%60);
        GLDimension *tDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: 1 + round(maxTime/sampleTime)  domainMin:0 length:round(maxTime/sampleTime)*sampleTime];
        tDim.name = @"time"; tDim.units = @"s";
        
        /************************************************************************************************/
        /*		Create the dynamical variables from the analytical solution								*/
        /************************************************************************************************/
        
        // We should check that we optimizations in places for these. They should be purely imaginary, and the multiplication and exponentiation should take advantage of that.
        GLFunction *iOmega = [[wave.eigenfrequencies swapComplex] makeHermitian];
        GLFunction *minusiOmega = [[[wave.eigenfrequencies swapComplex] negate] makeHermitian];
        
        NSArray * (^timeToUV) (GLScalar *) = ^( GLScalar *t ) {
            GLFunction *time_phase_plus = [[iOmega multiply: t] exponentiate];
            GLFunction *time_phase_minus = [[minusiOmega multiply: t] exponentiate];
            GLFunction *u = [[wave.Sprime transform: [[wave.u_plus multiply: time_phase_plus] plus: [wave.u_minus multiply: time_phase_minus]]] transformToBasis: @[@(kGLDeltaBasis), @(kGLDeltaBasis), @(kGLDeltaBasis)]];
            GLFunction *v = [[wave.Sprime transform: [[wave.v_plus multiply: time_phase_plus] plus: [wave.v_minus multiply: time_phase_minus]]] transformToBasis: @[@(kGLDeltaBasis), @(kGLDeltaBasis), @(kGLDeltaBasis)]];
            GLFunction *w = [[wave.S transform: [[wave.w_plus multiply: time_phase_plus] plus: [wave.w_minus multiply: time_phase_minus]]] transformToBasis: @[@(kGLDeltaBasis), @(kGLDeltaBasis), @(kGLDeltaBasis)]];
            GLFunction *rho = [[[[wave.S transform: [[wave.rho_plus multiply: time_phase_plus] plus: [wave.rho_minus multiply: time_phase_minus]]] transformToBasis: @[@(kGLDeltaBasis), @(kGLDeltaBasis), @(kGLDeltaBasis)]] times: wave.N2] plus: wave.rho];
            GLFunction *zeta = [[wave.S transform: [[wave.zeta_plus multiply: time_phase_plus] plus: [wave.zeta_minus multiply: time_phase_minus]]] transformToBasis: @[@(kGLDeltaBasis), @(kGLDeltaBasis), @(kGLDeltaBasis)]];
            
            return @[u,v,w,rho,zeta];
        };
        
        GLScalar *t = [GLScalar scalarWithValue: 0.0*2*M_PI/f0 forEquation: equation];
        NSArray *uv = timeToUV(t);
        GLFunction *u = uv[0];
        GLFunction *v = uv[1];
        GLFunction *w = uv[2];
        //GLFunction *rho = uv[3];
        GLFunction *zeta = uv[4];
        GLFunction *speed = [[[[u times: u] plus: [v times: v]] plus: [w times: w]] sqrt];
        GLFloat maxSpeed = [speed maxNow];
        NSLog(@"Initial maximum speed: %f", maxSpeed);
        
        /************************************************************************************************/
        /*		Let's also plop a float at a bunch of grid points.                                      */
        /************************************************************************************************/
        
        GLDimension *xFloatDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: ceil(width/horizontalFloatSpacingInMeters) domainMin: -width/2 length:width];
        xFloatDim.name = @"x-float";
        GLDimension *yFloatDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints:ceil(height/horizontalFloatSpacingInMeters)-1 domainMin: -height/2  length:height];
        yFloatDim.name = @"y-float";
        GLDimension *zFloatDim = [[GLDimension alloc] initWithPoints: @[ @(-38), @(-31.5), @(-25)]];
        zFloatDim.name = @"z-float";
        
        // For consistency, we order the float dimensions the same as the dynamical variable dimensions.
        NSArray *floatDimensions = @[zFloatDim, xFloatDim, yFloatDim];
        GLFunction *xFloat = [GLFunction functionOfRealTypeFromDimension: xFloatDim withDimensions: floatDimensions forEquation: equation];
        GLFunction *yFloat = [GLFunction functionOfRealTypeFromDimension: yFloatDim withDimensions: floatDimensions forEquation: equation];
        GLFunction *zFloat = [GLFunction functionOfRealTypeFromDimension: zFloatDim withDimensions: floatDimensions forEquation: equation];
        
        GLFunction *xIsopycnal = [GLFunction functionFromFunction: xFloat];
        GLFunction *yIsopycnal = [GLFunction functionFromFunction: yFloat];
        GLFunction *zIsopycnal = [GLFunction functionFromFunction: zFloat];
        GLFunction *isopycnalDeviation = [zeta interpolateAtPoints:@[zIsopycnal, xIsopycnal, yIsopycnal]];
        zIsopycnal = [zIsopycnal plus: isopycnalDeviation];
        
        GLFunction *xIsopycnalDiffusive = [GLFunction functionFromFunction: xIsopycnal];
        GLFunction *yIsopycnalDiffusive = [GLFunction functionFromFunction: yIsopycnal];
        GLFunction *zIsopycnalDiffusive = [GLFunction functionFromFunction: zIsopycnal];
        
        GLFunction *xFixedDepth = [GLFunction functionFromFunction: xFloat];
        GLFunction *yFixedDepth = [GLFunction functionFromFunction: yFloat];
        GLFunction *zFixedDepth = [GLFunction functionFromFunction: zFloat];
        
        GLFunction *xDrifter = [GLFunction functionFromFunction: xFloat];
        GLFunction *yDrifter = [GLFunction functionFromFunction: yFloat];
        GLFunction *zDrifter = [GLFunction functionFromFunction: zFloat];
        
        // Determine index ranges over which to average for the drogues.
        GLFloat drogueMin = -33;
        GLFloat drogueMax = -27;
        NSUInteger drogueMinIndex = NSNotFound;
        NSUInteger drogueMaxIndex = NSNotFound;
        for (NSUInteger iPoint=0; iPoint < zDim.nPoints; iPoint++) {
            if ([zDim valueAtIndex: iPoint] < drogueMin ) drogueMinIndex=iPoint;
            if ([zDim valueAtIndex: iPoint] < drogueMax ) drogueMaxIndex=iPoint;
        }
        NSRange drogueRange = NSMakeRange(drogueMaxIndex, drogueMaxIndex-drogueMinIndex+1);
        
        CGFloat cfl = 0.25;
        GLFloat cflTimeStep = cfl * xDim.sampleInterval / maxSpeed;
        GLFloat outputTimeStep = sampleTimeInMinutes*60;
        GLFloat timeStep = cflTimeStep > outputTimeStep ? outputTimeStep : outputTimeStep / ceil(outputTimeStep/cflTimeStep);
		
#warning Overriding time step!
		timeStep = 50;
		
        // Need this for the diffusive drifters
        GLFloat kappa = 5e-6; // m^2/s
        GLFloat norm = sqrt(timeStep*2*kappa);
        norm = sqrt(36./10.)*norm/timeStep; // the integrator multiplies by deltaT, so we account for that here.
        // RK4: dt/3 f(0) + dt/6 f(1) + dt/6 *f(4) + dt/3*f(3)
        // sqrt of ( (1/3)^2 + (1/6)^ + (1/6)^2 + (1/3)^2 )
        
        NSArray *floatArray=@[xIsopycnal, yIsopycnal, zIsopycnal, xIsopycnalDiffusive, yIsopycnalDiffusive, zIsopycnalDiffusive, xFixedDepth, yFixedDepth, xDrifter, yDrifter];
        GLAdaptiveRungeKuttaOperation *integrator = [GLAdaptiveRungeKuttaOperation rungeKutta4AdvanceY: floatArray stepSize: timeStep fFromTY:^(GLScalar *time, NSArray *yNew) {
            NSArray *uv = timeToUV(time);
            GLFunction *u2 = uv[0];
            GLFunction *v2 = uv[1];
            GLFunction *w2 = uv[2];
            
            GLFunction *xStep = [GLFunction functionWithNormallyDistributedValueWithDimensions: floatDimensions forEquation: wave.equation];
            GLFunction *yStep = [GLFunction functionWithNormallyDistributedValueWithDimensions: floatDimensions forEquation: wave.equation];
            GLFunction *zStep = [GLFunction functionWithNormallyDistributedValueWithDimensions: floatDimensions forEquation: wave.equation];
            xStep = [xStep times: @(norm)];
            yStep = [yStep times: @(norm)];
            zStep = [zStep times: @(norm)];
            
            GLSimpleInterpolationOperation *interpIso = [[GLSimpleInterpolationOperation alloc] initWithFirstOperand: @[u2, v2, w2] secondOperand: @[yNew[2], yNew[0], yNew[1]]];
            NSMutableArray *f = [interpIso.result mutableCopy];
            
            GLSimpleInterpolationOperation *interpIsoDiff = [[GLSimpleInterpolationOperation alloc] initWithFirstOperand: @[u2, v2, w2] secondOperand: @[yNew[5], yNew[3], yNew[4]]];
            NSArray *f2 = @[[interpIsoDiff.result[0] plus: xStep], [interpIsoDiff.result[1] plus: yStep], [interpIsoDiff.result[2] plus: zStep]];
            [f addObjectsFromArray: f2];
            
            GLSimpleInterpolationOperation *interpFixedDepth = [[GLSimpleInterpolationOperation alloc] initWithFirstOperand: @[u2, v2] secondOperand: @[zFixedDepth, yNew[6], yNew[7]]];
            NSArray *f3 = @[[interpFixedDepth.result[0] plus: xStep], [interpFixedDepth.result[1] plus: yStep]];
            [f addObjectsFromArray: f3];
            
            GLFunction *uMean = [uv[0] mean: 0 range: drogueRange];
            GLFunction *vMean = [uv[1] mean: 0 range: drogueRange];
            GLSimpleInterpolationOperation *interpDrifter = [[GLSimpleInterpolationOperation alloc] initWithFirstOperand: @[uMean, vMean] secondOperand: @[yNew[8], yNew[9]]];
            NSArray *f4 = @[[interpDrifter.result[0] plus: xStep], [interpDrifter.result[1] plus: yStep]];
            [f addObjectsFromArray: f4];
            
            return f;
        }];
        //integrator.absoluteTolerance = @[ @(1e0), @(1e0)];
        //integrator.relativeTolerance = @[ @(1e-6), @(1e-6), @(1e-6)];
        
        /************************************************************************************************/
        /*		Create a NetCDF file and mutable variables in order to record some of the time steps.	*/
        /************************************************************************************************/
        
//        NSString *outputFile = [[NSSearchPathForDirectoriesInDomains(NSDesktopDirectory, NSUserDomainMask, YES) objectAtIndex:0] stringByAppendingPathComponent: [NSString stringWithFormat: @"InternalWavesLatmix_%lu_%lu_%lu_GM_%.3f.nc", Nx, Ny, Nz_out,energyLevel]];
		NSString *outputFile = [NSString stringWithFormat: @"/Volumes/Data/InternalWavesLatmix_%lu_%lu_%lu_GM_%.3f.nc", Nx, Ny, Nz_out,energyLevel];
        GLNetCDFFile *netcdfFile = [[GLNetCDFFile alloc] initWithURL: [NSURL URLWithString: outputFile] forEquation: equation overwriteExisting: YES];
        
        [netcdfFile setGlobalAttribute: @(width) forKey: @"L_domain"];
        [netcdfFile setGlobalAttribute: @(latitude) forKey: @"latitude"];
        [netcdfFile setGlobalAttribute: @(f0) forKey: @"f0"];
        
        GLFunction *rhoScaled = [wave.rho scaleVariableBy: 1.0 withUnits: @"kg/m^3" dimensionsBy: 1.0 units: @"m"];
        rhoScaled.name = wave.rho.name;
        [netcdfFile addVariable: rhoScaled];
        
        GLFunction *n2Scaled = [wave.N2 scaleVariableBy: 1.0 withUnits: @"radians/s^2" dimensionsBy: 1.0 units: @"m"];
        n2Scaled.name = wave.N2.name;
        [netcdfFile addVariable: n2Scaled];
        
        integrator.shouldDisplayProgress = YES;
        [integrator integrateAlongDimension: tDim withTimeScale: 1.0 file: netcdfFile output: ^(GLScalar *t, NSArray *y) {
            NSMutableDictionary *scaledVariables = [NSMutableDictionary dictionary];
            NSArray *uv = timeToUV(t);
            
            scaledVariables[@"u"] = [uv[0] scaleVariableBy: 1.0 withUnits: @"m/s" dimensionsBy: 1.0 units: @"m"];
            scaledVariables[@"v"] = [uv[1] scaleVariableBy: 1.0 withUnits: @"m/s" dimensionsBy: 1.0 units: @"m"];
            scaledVariables[@"w"] = [uv[2] scaleVariableBy: 1.0 withUnits: @"m/s" dimensionsBy: 1.0 units: @"m"];
            scaledVariables[@"rho"] = [uv[3] scaleVariableBy: 1.0 withUnits: @"kg/m^3" dimensionsBy: 1.0 units: @"m"];
            scaledVariables[@"zeta"] = [uv[4] scaleVariableBy: 1.0 withUnits: @"m" dimensionsBy: 1.0 units: @"m"];
            scaledVariables[@"x-position"] = [y[0] scaleVariableBy: 1.0 withUnits: @"m" dimensionsBy: 1.0 units: @"float_id"];
            scaledVariables[@"y-position"] = [y[1] scaleVariableBy: 1.0 withUnits: @"m" dimensionsBy: 1.0 units: @"float_id"];
            scaledVariables[@"z-position"] = [y[2] scaleVariableBy: 1.0 withUnits: @"m" dimensionsBy: 1.0 units: @"float_id"];
            scaledVariables[@"x-position-diffusive"] = [y[3] scaleVariableBy: 1.0 withUnits: @"m" dimensionsBy: 1.0 units: @"float_id"];
            scaledVariables[@"y-position-diffusive"] = [y[4] scaleVariableBy: 1.0 withUnits: @"m" dimensionsBy: 1.0 units: @"float_id"];
            scaledVariables[@"z-position-diffusive"] = [y[5] scaleVariableBy: 1.0 withUnits: @"m" dimensionsBy: 1.0 units: @"float_id"];
            scaledVariables[@"x-position-fixed-depth"] = [y[6] scaleVariableBy: 1.0 withUnits: @"m" dimensionsBy: 1.0 units: @"float_id"];
            scaledVariables[@"y-position-fixed-depth"] = [y[7] scaleVariableBy: 1.0 withUnits: @"m" dimensionsBy: 1.0 units: @"float_id"];
            scaledVariables[@"z-position-fixed-depth"] = [zFixedDepth scaleVariableBy: 1.0 withUnits: @"m" dimensionsBy: 1.0 units: @"float_id"];
            scaledVariables[@"x-position-drifter"] = [y[8] scaleVariableBy: 1.0 withUnits: @"m" dimensionsBy: 1.0 units: @"float_id"];
            scaledVariables[@"y-position-drifter"] = [y[9] scaleVariableBy: 1.0 withUnits: @"m" dimensionsBy: 1.0 units: @"float_id"];
            scaledVariables[@"z-position-drifter"] = [zDrifter scaleVariableBy: 1.0 withUnits: @"m" dimensionsBy: 1.0 units: @"float_id"];
            
            return scaledVariables;
        }];
        
        [netcdfFile close];
    }
    return 0;
}
