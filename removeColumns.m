%% Assignment for AE4134: CFD I
% Created by:
% Zhi-li Liu 4146557
% Jasper van Wensveen 4142179
%
%removeColumns Removes the desired columns from a matrix
function [ strippedMatrix ] = removeColumns( matrix, columns )

% sort columns
sortedColumns = sort(columns);

% Get total amount of columns
nCols = size(matrix, 2);

% Set up vector with all indices of the columns
colIndicesTotal = 1:nCols;

% Determine columns to keep
strippedMatrixIndices = setdiff(colIndicesTotal, sortedColumns);

strippedMatrix = matrix(:,strippedMatrixIndices);
end

