%% Assignment for AE4134: CFD I
% Created by:
% Zhi-li Liu 4146557
% Jasper van Wensveen 4142179
%
%calculateAmountOfFluxes
function [ totalFluxes ] = calculateAmountOfFluxes( N )

% The amount of u-fluxes (horizontal fluxes) per row equals N+1. 
uFluxPerRow = N + 1;

% The amount of v-fluxes (vertical fluxes) per row is N.
vFluxPerRow = N;

% The amount of u-rows is N
amountURows = N;

% The amount of v-rows equals N+1
amountVRows = N + 1;

% Calculate total amount of fluxes
totalFluxes = uFluxPerRow*amountURows + vFluxPerRow*amountVRows;

end

