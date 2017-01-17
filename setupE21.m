%% Assignment for AE4134: CFD I
% Created by:
% Zhi-li Liu 4146557
% Jasper van Wensveen 4142179
%
%setupE21 Setup the incidence matrix E21, which is an extended matrix.
%(Added extra bounds for the inner-grid to close it)
function [ E21 , totalFluxes] = setupE21( N )


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

% Due to addition of the closing boundaries, the amount of 2-cells is
% increased by 1 in each direction. (Actually 2, but N was defined for the
% outer-grid, which has 1 cell more than the inner-grid)
newN  = N + 1;

% The incidence marix would look something like this per row [.. +1 .. -1
% +1 ... -1]

% Create the +1's and -1's
firstDiag = ones(newN^2, 1);
secondDiag = -ones(newN^2, 1);
thirdDiag = ones(newN^2, 1);
fourthDiag = -ones(newN^2, 1);

% Assemble into s vector
s = [firstDiag; secondDiag; thirdDiag; fourthDiag];

% The i indices indicate which rows are used. Here we use all rows 4 times.
% Total amount of rows equals newN^2
i = [1:newN^2, 1:newN^2, 1:newN^2, 1:newN^2]';

% Determine the columns for each of the 'diagonals' (might not be real
% diagonals, rather partially diagonals)

% The first J is very simple. Just starting from 1 to the amount of cells.
% This is because the horizontal fluxes (u) are dealt with first in the u
% vector.
firstJ = (1:newN^2);

% The second J is the vertical component on the left side of a 2-cell. The
% first element starts at totalFluxes/2 +1, then has a few elements
% sequentially after it, jumps one and this pattern is continued.

% Get all vertical flux (v) indices
verticalFluxIndices = (totalFluxes/2 +1):totalFluxes;

% Remove the last flux in each row (vFluxPerRow)
removeElementIndicator = (mod(verticalFluxIndices, vFluxPerRow) == 0);
keepElementIndicator = not(removeElementIndicator);

secondJ = verticalFluxIndices(keepElementIndicator);

% The third J is just one shifted from the second one
thirdJ = secondJ + 1;

% The fourth J is just one row (uFluxPerRow) shifted from firstJ
fourthJ = firstJ + uFluxPerRow;

% assemble j matrix
j = [firstJ, secondJ, thirdJ, fourthJ]';

% create sparse matrix
E21 = sparse(i,j,s);
end

