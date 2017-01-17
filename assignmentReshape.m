%% Assignment for AE4134: CFD I
% Created by:
% Zhi-li Liu 4146557
% Jasper van Wensveen 4142179
%
%assignmentReshape Reshapes a vector into a matrix, using the convention
%used in the CFD I assignment. That is: from left to right, from bottom to
%top. So a vector [1, 2 ..... ,25] reshaped into a 5x5 matrix, should have
%at the top row [21, 22, 23, 24, 25]
%
% The input parameter vector indicates the input vector, this is the vector
% to be reshaped.
% n indicates the amount of rows and m indicates the amount of columns
%
% A function is written, as this provides more flexibility when the reshape
% changes.

function [ output ] = assignmentReshape( vector, n, m )

output = reshape(vector, n, m)';
end

