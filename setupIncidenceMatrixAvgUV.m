%% Assignment for AE4134: CFD I
% Created by:
% Zhi-li Liu 4146557
% Jasper van Wensveen 4142179
%
%setupIncidenceMatrixAvgUV

function [ U, V ] = setupIncidenceMatrixAvgUV( N )

tE21 = setupTE21(N, calculateAmountOfFluxes(N));
tE21 = full(tE21);

U = tE21(:,1:end/2);
V = tE21(:,end/2+1:end);

% Now turn them into 0.5
U = abs(U)*0.5;
V = abs(V)*0.5;
end

