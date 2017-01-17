%% Assignment for AE4134: CFD I
% Created by:
% Zhi-li Liu 4146557
% Jasper van Wensveen 4142179
%
% boundaryUIndices This grabs the outer oriented matrix and determines which
%of these indices are boundary ones.
% The numbering scheme was bottom to top, left to right, first u then v
function [ boundaryIndices ] = boundaryUIndices( N )


amountOfFluxes = calculateAmountOfFluxes( N );

% Fluxes per direction
amountOfFluxesOneDir = amountOfFluxes/2;

% The amount of u fluxes per row is
uFluxPerRow = amountOfFluxesOneDir/N;

% The amount of v fluxes per column
vFluxPerCol = N;

% Left boundary u indices
leftIndices = (0:N-1).*uFluxPerRow+1;

% Right boundary u indices
rightIndices = leftIndices + uFluxPerRow -1;

% Bottom boundary u (v) indices
bottomIndices = (amountOfFluxesOneDir+1:amountOfFluxesOneDir+vFluxPerCol);

% Top boundary u (v) indices
topIndices = (amountOfFluxes - vFluxPerCol + 1 : amountOfFluxes);

% Assemble the final output vector.
boundaryIndices = [leftIndices, rightIndices, bottomIndices, topIndices]';
end

