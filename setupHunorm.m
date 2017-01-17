%% Assignment for AE4134: CFD I
% Created by:
% Zhi-li Liu 4146557
% Jasper van Wensveen 4142179
%
function [Hu_norm] = setupHunorm(N, Hu_norm, th, h)

% Get boundaryIndices
boundaryIndices = boundaryUIndices( N );

% Create boundary prescribed H1t1 in the sequence [left, right, bottom, top]                
prescribedH1t1 = [  h(1)./th';...
                    h(end)./th';...
                    h(1)./th';...
                    h(end)./th']; 
                
% calculate the Hu_norm
Hu_norm(boundaryIndices) = prescribedH1t1;
