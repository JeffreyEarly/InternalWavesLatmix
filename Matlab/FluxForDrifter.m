% Given a vector y where size(y)=[N 8] such that we have N floats in three
% different experiments. They are ordered [x_float y_float z_float x_diff
% y_diff z_diff x_drifter y_drifter], and should therefore return [u v w u
% v w u v].
function flux = FluxForDrifter(t,y,z_drifter,wavemodel,method)    
    [u,v] = wavemodel.VelocityAtTimePosition(t,y(:,1),y(:,2),z_drifter, method);
    flux = cat(2,u,v);
end