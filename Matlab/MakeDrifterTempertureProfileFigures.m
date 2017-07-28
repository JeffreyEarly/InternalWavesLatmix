file = '/Volumes/OceanTransfer/FloatsWithTemperatureProfileExperiment_2017-07-26T221153_256x64x32.nc';

t = ncread(file, 't');
xCoord = ncread(file, 'x');
yCoord = ncread(file, 'y');
zCoord = ncread(file, 'z');
N2 = ncread(file, 'N2');
Nx = length(xCoord);
Ny = length(yCoord);
Nz = length(zCoord);

x = ncread(file, 'x-position-drifter')';
y = ncread(file, 'y-position-drifter')';
z = ncread(file, 'z-position-drifter')';

% [x_com, y_com, q, r] = CenterOfMass( x, y );
% 
% figure, plot(q,r)

iDrifter = 1;
rhoAtDrifter = zeros(length(t),Nz);
[X,Y,Z] = ndgrid(xCoord,yCoord,zCoord);
for iTime = 1:length(t)
    rho3d = squeeze(ncread(file, 'rho', [1 1 1 iTime], [Nx Ny Nz 1], [1 1 1 1]));
    [xq,yq,zq] = ndgrid(x(iTime,iDrifter),y(iTime,iDrifter),zCoord);
    rhoAtDrifter(iTime,:) = interpn(X,Y,Z,rho3d,xq,yq,zq);
end

[Tt,Zt] = ndgrid(t,zCoord);
depth = single(linspace(-35,-10,100));
[T,D] = ndgrid(t,depth);
rhoAtDepth = interpn(Tt,Zt,rhoAtDrifter,T,D);

% These numbers are fairly random
rho0 = 1027;
T0 = 10;
alpha = 0.15; % expansion coefficient for the equation of state kg/(m^3 C)
temperatureAtDepth = T0 - (rhoAtDepth - rho0)/alpha;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Draw the figure
%

% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
run('LoadFigureDefaults')

figure('Units', 'points', 'Position', [50 50 1000 300])
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Draw buoyancy frequency vs depth
%
ax1 = subplot(1,5,1);
plot(sqrt(N2)*3600/(2*pi),zCoord, 'Color', 'black', 'LineWidth', 2)
xlabel('cycles/hr', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('z (meters)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylim([min(depth) max(depth)])

ax2 = subplot(1,5,[2 3 4 5]);
pcolor(t/86400,depth,temperatureAtDepth'), shading flat, colormap(jet), hold on
C = contourc(double(t),double(depth(depth<-20)),double(temperatureAtDepth(:,depth<-20))',[19 19]);
i = 1;
while i < size(C,2)
    j = i+C(2,i); % end
    xc = C(1,(i+1):j);
    yc = C(2,(i+1):j);
    plot(xc/86400,yc, 'Color', 'black', 'LineWidth', 2);
    i = j+1;
end
xlabel('time (days)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
set(gca,'YTickLabel',[])

c = colorbar('eastoutside');
c.Label.String = 'T (°C)';
caxis([17 23])

ax2.Position = [ax2.Position(1)-0.03 ax2.Position(2) ax2.Position(3) ax2.Position(4)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Finally, draw the schematic drifter on the first plot, now that sizes
% 	are fixed.
%
axes(ax1)

drogueDepth = 27;
drogueHeight = 6;
drifterSigma = 20; % x offset in units of sigma (non physical)
drifterWidthSigma = 8;

r = ax1.PlotBoxAspectRatio(1) * ax1.DataAspectRatio(2);
tetherWidthSigma = 0.05*drifterWidthSigma;
drogueWidthSigma = 0.4*drifterWidthSigma;

rectangle('Position',[drifterSigma-0.5*tetherWidthSigma -drogueDepth-1 tetherWidthSigma drogueDepth+2], 'FaceColor', 0.5*ones(3,1))
ball = rectangle('Position',[drifterSigma-0.5*drifterWidthSigma -1 drifterWidthSigma r*drifterWidthSigma],'Curvature',[1 1], 'FaceColor', 'yellow', 'Clipping', 'off');
% Little hack to get clipping to work correctly.
ball.Clipping = 'on';
pause(.5)
ball.Clipping = 'off';
pause(.5)
drogue = rectangle('Position',[drifterSigma-0.5*drogueWidthSigma -drogueDepth-drogueHeight drogueWidthSigma drogueHeight],'Curvature',[0.1 0.1], 'FaceColor', 'yellow', 'FaceColor', 'black');

print('-dpng', '-r300', 'DrifterTemperatureProfile.png')