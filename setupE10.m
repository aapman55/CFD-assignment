%% Assignment for AE4134: CFD I
% Created by:
% Zhi-li Liu 4146557
% Jasper van Wensveen 4142179
%
%setupE10
function [ E10 ] = setupE10( N )


% Get the amount of circulations (on the inner grid this is called the
% circulation)
amountOfCirc = calculateAmountOfFluxes(N);

% Get rid of the boundaries. Each side of the square has N boundary fluxes.
amountOfCircwithoutBoundaries = amountOfCirc - 4 * N;

% Calculate the amount of sinks
amountOfSinks = N^2;

% Get all the indices of the sinks
sinkIndices = 1:amountOfSinks;

% For the horizontal circulation, there are 2 diagonals (sort of, with jump
%s between rows). For these 2 diagonals we are going to determine the
%corresponding sink.

% Start with the sink on the left side of the horizontal circulation
removeIndices = (mod(sinkIndices, N) == 0);
keepIndices = not(removeIndices);

% The keepIndices is a boolean operator and indicates which elements in the
% sinkIndices are going to be used.
sinkLeftHorizontalIndices = sinkIndices(keepIndices);

% The sink on the right side of the horizontal circulation is just one
% index number higher than the left.
sinkRightHorizontalIndices = sinkLeftHorizontalIndices + 1;

% Now let's look at the vertical circulations. For the vertical case we
% also get 2 kind of diagonals. One at the bottom of the circulation and
% one at the top of the circulation.

% First we do the bottom ones.
% Create a vector containing all the indices of the vertical circulations.
% We look at the vertical circulations only, this means that we are
% starting the count at 1. By doing so, the formula to get from the
% circulation number to the desired sink number is simplified.
verticalCircNumbers = 1:amountOfCircwithoutBoundaries/2;

% The bottom sink correspond just to the vertical circulation number.
sinkBottomVerticalIndices = verticalCircNumbers;

% The top one differs one row. Meaning that there is an increase in the
% index number equal to the amount of sinks on one row. This value is equal
% to N.
sinkTopVerticalIndices = sinkBottomVerticalIndices + N;

% At this point, we have all the column indices for the E10 matrix. The row
% indices just run from the first row to halfway or from halfway to the end, 
% depending on the type of circulation (horizontal or vertical). The
% horizontal circulations will occupy the first half of the E10 matrix.
% The order is : left sink, right sink, bottom sink, top sink.
i = [   1:amountOfCircwithoutBoundaries/2;...
        1:amountOfCircwithoutBoundaries/2;...
        (amountOfCircwithoutBoundaries/2 + 1):amountOfCircwithoutBoundaries;...
        (amountOfCircwithoutBoundaries/2 + 1):amountOfCircwithoutBoundaries];

% Assemble the j vector for the column numbers
j = [sinkLeftHorizontalIndices; sinkRightHorizontalIndices; sinkBottomVerticalIndices; sinkTopVerticalIndices];

% Now give the location (i,j) in matrix E10 a value. The left and bottom
% will be -1 and the right and top will be 1. This is due to the sign
% convention that there are sinks. So moving into the sink will be
% considered as positive.
s = [   -ones(size(sinkLeftHorizontalIndices));...
        ones(size(sinkRightHorizontalIndices));...
        -ones(size(sinkBottomVerticalIndices));...
        ones(size(sinkTopVerticalIndices))];
    
% Assemble the sparse matrix E10;
E10 = sparse(i, j, s);
end

