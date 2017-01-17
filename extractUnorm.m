%% Assignment for AE4134: CFD I
% Created by:
% Zhi-li Liu 4146557
% Jasper van Wensveen 4142179
%
%extractUnorm Using the set boundary conditions, calculate the u_norm
%vector. 
function [ u_norm ] = extractUnorm( tE21, U_wall_left, U_wall_right, V_wall_bot, V_wall_top, th )


% Extract N2 from tE21
N2 = size(tE21,1);

% Get N
N = sqrt(N2);

% Get boundaryIndices
boundaryIndices = boundaryUIndices(N);

% Create boundary prescribed u in the sequence [left, right, bottom, top]
prescribedU = [ ones(N,1)*U_wall_left.*th';...
                ones(N,1)*U_wall_right.*th';...
                ones(N,1)*V_wall_bot.*th';...
                ones(N,1)*V_wall_top.*th'];

% calculate the u_norm
u_norm = tE21(:,boundaryIndices)*prescribedU;

end

