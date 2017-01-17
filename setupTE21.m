%% Assignment for AE4134: CFD I
% Created by:
% Zhi-li Liu 4146557
% Jasper van Wensveen 4142179
%
%setupTE21 Setup the incidence matrix that goes from the outer oriented
%grid to the inner oriented grid. Relating the fluxes to the volumes.
%
% Numbering scheme:
% From bottom to top, left to right, first u then v.
function [ tE21 ] = setupTE21( N, nFluxes )

% The amount of fluxes in one direction
nFluxOneDir = nFluxes/2;

% The amount of cells in each direction
N2 = N^2;

% amount of u fluxes per row
uFluxPerRow = N + 1;

% The incidence matrix will contain 4 diagonals. 2 of which are -1 and the
% other 2 are 1. The correct sequence per row is [-1 1 ... -1 ... 1]
firstDiag = -ones(N2,1);
secondDiag = ones(N2,1);
thirdDiag = -ones(N2,1);
fourthDiag = ones(N2,1);

s = [firstDiag; secondDiag; thirdDiag; fourthDiag];

% Now determine the i indices (i index indicates the row numbers) These
% diagonals span all rows.
oneRow = (1:N2)';

i = [oneRow; oneRow; oneRow; oneRow];

% Now determine the j indices (column numbers) Starting points are 1 2 N2+1
% and N2+1+N

% The first and second J element have a jump in them

% Get all the indices of the horizontal elements
horizontalFluxIndices = 1:nFluxOneDir;

% Remove the last flux in each row (vFluxPerRow)
removeElementIndicator = (mod(horizontalFluxIndices, uFluxPerRow) == 0);
keepElementIndicator = not(removeElementIndicator);

firstJ = horizontalFluxIndices(keepElementIndicator)';

j = [firstJ; firstJ+1; oneRow+nFluxOneDir; oneRow+nFluxOneDir+N];

tE21 = sparse(i,j,s);

end

