%% Assignment for AE4134: CFD I
% Created by:
% Zhi-li Liu 4146557
% Jasper van Wensveen 4142179
%
% setupConvectionIncidenceMatrices
%
% The convective term looks something like (as written in assignment):
%
% - h2/4*h2 *(v1,2 + v2,2)* xi2,2 - h2/4*h3 * (v1,3 + v2,3) * xi2,3
%
% However, in the code it is first upper, than lower, and first left than right side.
%
% Matrices:
% Matrix to determine the vector corresponding to h2 (denominator) in this example    [1]
% Matrix to determine the vector corresponding to h3 (denominator) in this example    [2]
% Matrix to determine (v1,2 + v2,2) in this example                     [3]
% Matrix to determine xi2,2 in this example                             [4]
% Matrix to determene (v1,3 + v2,3) in this example                     [5]
% Matrix to determine xi2,3 in this example                             [6]
% Matrix to determine h2 in the nominator                               [7]
%
% Matrix 7 indicates the length in the direction of the unknown flux
%
% Please note that this example consist of just one row. In reality there
% will be multiple rows
%
% Basically you will get:
%
% - ([7]*{h})./(4*[1]*{h}) .* ([3]*{u}) .* ([4]*{xi}) -
% ([7]*{h})./(4*[2]*{h}) .* ([5]*{u}) .* ([6]*{xi})
%
% Where {u} indicates the u vector including all fluxes (including the
% extended ones). Each iteration the unknown u's in this {u} vector needs to
% be updated
%
% {h} indicates the vector containing all the lengths
% {xi} containing all the xis
function [ M1, M2, M3, M4, M5, M6, M7 ] = setupConvectionIncidenceMatrices( N )


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

% Determine the amount of h (length of inner grid elements)
amountH = uFluxPerRow;

% Determine amount xi per row
xiPerRow = amountH;

% Determine amount of xis
amountXi = xiPerRow^2;

% Get indices for the extended boundaries
extendedBoundaryIndices = extendedBoundaryUIndices(N);

% Get the indices of the unknown u flows (u and v)
unknownFlowIndices = setdiff(1:totalFluxes, extendedBoundaryIndices);

% Initiate Matrices
% M1 = zeros(length(unknownFlowIndices), uFluxPerRow);
% M2 = zeros(length(unknownFlowIndices), uFluxPerRow);
% 
% M3 = zeros(length(unknownFlowIndices), totalFluxes);
% M4 = zeros(length(unknownFlowIndices), amountXi);
% 
% M5 = zeros(length(unknownFlowIndices), totalFluxes);
% M6 = zeros(length(unknownFlowIndices), amountXi);
% 
% M7 = zeros(length(unknownFlowIndices), uFluxPerRow);

% Initialise indices for the matrices
indexM1 = zeros(length(unknownFlowIndices),1);
indexM2 = zeros(length(unknownFlowIndices),1);
indexM3_1 = zeros(length(unknownFlowIndices),1);
indexM3_2 = zeros(length(unknownFlowIndices),1);
indexM4 = zeros(length(unknownFlowIndices),1);
indexM5_1 = zeros(length(unknownFlowIndices),1);
indexM5_2 = zeros(length(unknownFlowIndices),1);
indexM6 = zeros(length(unknownFlowIndices),1);
indexM7 = zeros(length(unknownFlowIndices),1);

% Initialise values for the matrices
valM1 = zeros(length(unknownFlowIndices),1);
valM2 = zeros(length(unknownFlowIndices),1);
valM3_1 = zeros(length(unknownFlowIndices),1);
valM3_2 = zeros(length(unknownFlowIndices),1);
valM4 = zeros(length(unknownFlowIndices),1);
valM5_1 = zeros(length(unknownFlowIndices),1);
valM5_2 = zeros(length(unknownFlowIndices),1);
valM6 = zeros(length(unknownFlowIndices),1);
valM7 = zeros(length(unknownFlowIndices),1);

% Iterate over the unknown fluxes
for i = 1:length(unknownFlowIndices)
    uNumber = unknownFlowIndices(i);
    % The u fluxes
    if (i <= length(unknownFlowIndices)/2)
        URow = ceil(uNumber/uFluxPerRow);
        UCol = mod(uNumber, uFluxPerRow);
        
        % Index for Matrix 1
        indexM1(i) = URow;

        % Index for Matrix 2
        indexM2(i) = indexM1(i) - 1;
        
        % Indices for Matrix 3
        indexM3_1(i) = 0.5*totalFluxes + vFluxPerRow*(URow - 1) + UCol;
        indexM3_2(i) = indexM3_1(i) + 1;
        
        % Index for Matrix 4
        indexM4(i) = (URow - 1)*xiPerRow + UCol;
        
        % Indices for Matrix 5
        indexM5_1(i) = indexM3_1(i) - vFluxPerRow;
        indexM5_2(i) = indexM3_2(i) - vFluxPerRow;
        
        % Index for Matrix 6
        indexM6(i) = indexM4(i) - xiPerRow;
        
        % Index for Matrix 7
        indexM7(i) = UCol;
       
    % The v fluxes
    else
        VRow = ceil( (uNumber-0.5*totalFluxes) / vFluxPerRow);
        VCol = mod( (uNumber-0.5*totalFluxes) , vFluxPerRow );
        
        % Index for Matrix 1
        indexM1(i) = VCol - 1;

        % Index for Matrix 2
        indexM2(i) = indexM1(i) + 1;
        
        % Indices for Matrix 3
        indexM3_1(i) = uFluxPerRow*(VRow - 1)+ VCol - 1;
        indexM3_2(i) = indexM3_1(i) + uFluxPerRow;
    
        % Indices for Matrix 4
        indexM4(i) = (VRow - 1)*xiPerRow + VCol - 1;
        
        % Indices for Matrix 5
        indexM5_1(i) = indexM3_1(i) + 1;
        indexM5_2(i) = indexM3_2(i) + 1;
        
        % Indices for Matrix 6
        indexM6(i) = indexM4(i) + 1;
        
        % Indices for Matrix 7
        indexM7(i) = VRow;

    end
    
    % Set the sign of the xi's. This is done as the u's have negative signs
    % in front of the equation and v's positive signs. To follow the
    % notation in the assignment (i.e. the example given in the assignment
    % for u), a negative sign will be given to the v xi's. In the main
    % there will be a minus sign in front of each convection part, making
    % the v's positive again.
    xiSign = 1;
    
    if ((i > length(unknownFlowIndices)/2))
        xiSign = -1;
    end
    
    % Set indices Matrix M1
    valM1(i) = 1;
    
    % Set indices Matrix M2
    valM2(i) = 1;
    
    % Set incidence matrix M3
    valM3_1(i) = 1;
    valM3_2(i) = 1;
    
    % Set indices Matrix M4
    valM4(i) = 1 * xiSign;

    % Set indices matrix M5
    valM5_1(i) = 1;
    valM5_2(i) = 1; 
    
    % Set indices Matrix M6
    valM6(i) = 1 * xiSign;
    
    %Set indices Matrix M7
    valM7(i) = 1;
    
end
counter = (1:length(unknownFlowIndices))';

% Matrices M3, M4, M5, M6 and M7 end with zero columns. Meaning that in the
% sparse matrix they will not be there. Resulting in a erroneous size of
% the matrix. To overcome this the bottom right entry will be defined as 0.

% Construct the sparse matrices
M1 = sparse(counter, indexM1, valM1);
M2 = sparse(counter, indexM2, valM2);
M3 = sparse([counter;counter;counter(end)], [indexM3_1;indexM3_2;totalFluxes], [valM3_1;valM3_2;0]);
M4 = sparse([counter;counter(end)], [indexM4; amountXi], [valM4;0]);
M5 = sparse([counter;counter;counter(end)], [indexM5_1;indexM5_2;totalFluxes], [valM5_1;valM5_2;0]);
M6 = sparse([counter;counter(end)], [indexM6; amountXi], [valM6;0]);
M7 = sparse([counter; counter(end)], [indexM7; uFluxPerRow], [valM7;0]);

end

