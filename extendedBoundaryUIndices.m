%% Assignment for AE4134: CFD I
% Created by:
% Zhi-li Liu 4146557
% Jasper van Wensveen 4142179
%
%extendedBoundaryUIndices For the extended matrix, extract all the indices
%of the boundaries (both tangential as normal)
%
% The numbering scheme was bottom to top, left to right, first u then v
%
% N is the REAL N. Thus coorresponding to the outer-grid (non-extended)
function [ allIndices ] = extendedBoundaryUIndices( N )


% The amount of u-fluxes (horizontal fluxes) per row remains the same and
% equals N+1. 
uFluxPerRow = N + 1;

% The amount of v-fluxes (vertical fluxes) per row increases by 2. Since on
% both sides a closing flux is inserted. Normally the amount is N, now it
% is N+2.
vFluxPerRow = N + 2;

% Now determine the amount of u-rows present. Due to the addition of a
% closing top and bottom u flux, the amount of rows for the u component
% will be 2 more than normal. Normally the amount of u-rows is N. Now it is
% N+2.
amountURows = N + 2;

% Determine the amount of rows for the v-fluxes. This amount does not
% change and remains N+1
amountVRows = N + 1;

% Calculate total amount of fluxes
totalFluxes = uFluxPerRow*amountURows+vFluxPerRow*amountVRows;

% First perform analysis on extended boundaries (tangential fluxes)

% Bottom tangential Indices
tangentialBottom = 1:uFluxPerRow;

% Top tangential indices
tangentialTop = (totalFluxes/2 - uFluxPerRow + 1):totalFluxes/2;

% Left tangential indices
tangentialLeft = (0:amountVRows-1)*vFluxPerRow+totalFluxes/2+1;

% Right tangential indices
tangentialRight = (1:amountVRows)*vFluxPerRow+totalFluxes/2;

% Now move on to the normal components

% Left normal component
normalLeft = (uFluxPerRow + 1) + (0:amountURows - 2-1)*uFluxPerRow;

% Right normal component
normalRight = normalLeft + (uFluxPerRow - 1);

% bottom normal component
normalBottom = totalFluxes/2+1+(1:vFluxPerRow - 2);

% Top normal components
normalTop = totalFluxes - vFluxPerRow + 2:totalFluxes -1;

% Create output vector
allIndices = [  tangentialBottom, tangentialTop, tangentialLeft, tangentialRight, ...
                normalLeft, normalRight, normalBottom, normalTop];

% sort output vector
allIndices = sort(allIndices)';
end

