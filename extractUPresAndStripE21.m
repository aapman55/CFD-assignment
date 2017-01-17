%% Assignment for AE4134: CFD I
% Created by:
% Zhi-li Liu 4146557
% Jasper van Wensveen 4142179
%
%extractUPresAndStripE21 Use known boundary conditions to create a vector
%with known values (UPres). Then at the same time remove the corresponding
%columns in E21.
function [ UPres, E21 ] = extractUPresAndStripE21( E21, U_wall_top, U_wall_bot, U_wall_left, U_wall_right, ...
                                            V_wall_top, V_wall_bot, V_wall_left, V_wall_right , h)

% Extract current Extended N^2
ExtendedN2 = size(E21, 1);

% Determine Extended N
ExtendedN = sqrt(ExtendedN2);

% Real N
N = ExtendedN - 1;

% Get vector for the prescribed u velocities at the boundaries
% BEWARE: For this function the REAL N is used, not the extended!
uPresVal = setupupres(N, U_wall_top, U_wall_bot, U_wall_left, U_wall_right, ...
                                            V_wall_top, V_wall_bot, V_wall_left, V_wall_right, h);
                                        
% Get indices for the E21 index that correspond to the prescribed u's
prescribedIndices  = extendedBoundaryUIndices( N );

% Get the UPres, which is an output of this function
UPres = E21(:,prescribedIndices)*uPresVal;

% Remove the columns in E21 that has been prescribed
E21 = removeColumns(E21, prescribedIndices);

end

