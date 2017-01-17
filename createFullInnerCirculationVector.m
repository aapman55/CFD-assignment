%% Assignment for AE4134: CFD I
% Created by:
% Zhi-li Liu 4146557
% Jasper van Wensveen 4142179
%
%createFullInnerCirculationVector Creates a vector with the right amount of circulation and
%also fills in the prescribed boundary circulations. The output will be the
%vector with the indices that still needs to be filled in.
% The assumption made is that the amount of boundary circulations are the same
% for left, right, bottom and top. And the assumed numbering is from left
% to right bottom to top first u then v.
function [ fluxVector, unknownUIndices ] = createFullInnerCirculationVector( N, U_wall_left, U_wall_right, V_wall_bot, V_wall_top  )

% Calculate the amount of fluxes
amountOfFluxes = calculateAmountOfFluxes(N);

% Get the indices for the boundaries
boundaryIndices = boundaryUIndices(N);

% Amount of fluxes on the boundary
amountOfBoundaryFluxes = length(boundaryIndices);

% Amount of boundary fluxes per side
amountOfBoundaryFluxesPerSide = amountOfBoundaryFluxes/4;

% Create flux vector
fluxVector = zeros(amountOfFluxes, 1);

% Create dummy vector
dummy = ones(amountOfBoundaryFluxesPerSide, 1);

% Fill in the boundary
fluxVector(boundaryIndices) = [U_wall_left*dummy, U_wall_right*dummy, V_wall_bot*dummy, V_wall_top*dummy];

% Create the vector with indices to fill in
unknownUIndices = setdiff(1:amountOfFluxes, boundaryIndices);
end

