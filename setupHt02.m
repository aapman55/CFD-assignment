%% Assignment for AE4134: CFD I
% Created by:
% Zhi-li Liu 4146557
% Jasper van Wensveen 4142179
%
% Mapping the inner-oriented 2-cochains to the outer-oriented 0-cochains means dividing by the
% area of the square.
function [Ht02] = setupHt02(N, h)

% Create an array of all the inner-oriented 2-cochains dx
tdx = repmat(h,1,N+1);

% Create an array of all the inner-oriented 2-cochains dy
tdy = reshape(repmat(h,N+1,1),1,[]);

% Set-up Ht02 by dividing 1 by the multiplication of the dx and dy (areas)
Ht02 = sparse(diag(1./(tdx.*tdy)));
