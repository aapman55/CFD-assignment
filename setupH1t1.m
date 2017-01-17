function [H1t1] = setupH1t1(N, th, h)
% tu/tdeltaY = u/deltaX  => u = tu*deltaX/tdeltaY;
% where tdeltaY is the length of the outer grid segment and deltaX the length
% of the inner grid segment

% Set-up the array containing the flux (outer-grid) area lengths for each velocity
outer = [reshape(repmat(th,N+1,1),1,[]),repmat(th,1,N+1)];

% Set-up the array containing the vorticity (inner-grid) area lengths for each velocity
inner = [repmat(h,1,N),reshape(repmat(h,N,1),1,[])];

% Construct the diagonal Hodge matrix
% H1t1 = sparse(diag(inner./outer)); % vroegere beun
i = 1:length(inner);
H1t1 = sparse(i,i,inner./outer);
