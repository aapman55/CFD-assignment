%% Assignment for AE4134: CFD I
% Created by:
% Zhi-li Liu 4146557
% Jasper van Wensveen 4142179

close all
clear
clc

% This system that you need to solve will be singular. MATLAB gives you a
% warning at each time step. To switch of this warning, remove the comment
% in the next line

warning off

% This file contains the skeleton of the program which will solve the lid
% driven cavity problem on a unit square. The parts that have to be
% supplemented are described in the assignment.
%
% The pieces that need to be filled in are indicated by dots: ............
%

%
% When running the code, determine a suitable time step. A too small time
% step will make the calculation very long, while for a too large time step
% the solution will blow up due to numerical instability.
%
tic;                   % Start timing
Re = 1000;             % Reynolds number
N = 32;                 % Number of volumes in the x- and y-direction
Delta = 1/N;           % uniforme ruimte stap in de afbeelding
% The time step is calculated later based on the grid size.
t = 0.0001;            % time (0.0001)
tol = 1e-6;            % tol determines when steady state is reached and the program terminates
turbo = 5;             % A value between 1 and 5 to increase speed



% wall velocities
U_wall_top = -1;
U_wall_bot = 0;
U_wall_left = 0;
U_wall_right = 0;
V_wall_top = 0;
V_wall_bot = 0;
V_wall_left = 0;
V_wall_right = 0;

%%
%   Generation of a non-uniform mesh
%

%
%   tx are the coordinates of the nodal points on the outer-oriented grid
%
tx = zeros(1,N+1);
for i=1:N+1
    txi = (i-1)*Delta;
    tx(i) = 0.5*(1 - cos(pi*txi));
end

% Local mesh size on outer oriented grid
th = tx(2:N+1) - tx(1:N);

%
% x are the coordinates of the nodal points on the inner-oriented grid (including
% endpoints 0 and 1)
% h contains the edge lengths on the inner-oriented grid
%
x = 0.5*(tx(1:N) + tx(2:N+1));
x = [0 x 1];

h = x(2:N+2) - x(1:N+1);


%%%%%%%%% Recalculate dt based on minimum criteria
dt = turbo*min(0.5*Re*min(h)^2,min(h));


%   Initial condition u = v = 0
%
%   Both u and v will be stored in one big vector called 'u'
%
%   The vector u only contains the true unknowns, not the velocities that
%   are prescribed by the boundary conditions
%
%   The vector u contains the *inner-oriented* fluxes as unknowns
%
%u = ones(2*N*(N-1),1); % Originally given in the assignment, but does not correspond to
%comments
u = zeros(2*N*(N-1),1);



%%
% Set up the Incindence matrix 'tE21' which connects the fluxes to the
% volumes. Use the orientation described in the assignment.
%
tE21 = setupTE21(N, calculateAmountOfFluxes(N));
%ftE21 = full(tE21);

%
% Inserting boundary conditions for normal velocity components
%                                                        
u_norm = extractUnorm( tE21, U_wall_left, U_wall_right, V_wall_bot, V_wall_top , th);% outer/primal
%fu_norm = full(u_norm);

%
% Remove columns associated with prescribed normal velocities from
% Incidenc,e matrix 'tE21'
%
tE21 = removeColumns(tE21, boundaryUIndices(N));
%ftE21 = full(tE21)

% ==============================================
% Check the tE21 matrix using E10
% ==============================================
E10 = setupE10(N);
if (prod(prod(E10 == -tE21')) ~= 1)
    error('Please check your tE21 matrix!');
end

%%
% Setting up simple Hodge matrix which converts the fluxes on the
% outer-oriented grid to circulation on the inner-oriented grid. Assume
% that the fluxes and circulation are constant over each 1-cell. This will
% give a diagonal Hodge matrix. Call this Hodge matrix 'H1t1'
%
H1t1 = setupH1t1(N, th, h);

%
% Pre-allocate Hu_norm vector (was already present in skeleton)
%
Hu_norm = zeros(2*N*(N+1),1);

%
% Hu_norm is the vector which will contain the Hodge of the prescribed
% normal fluxes. Calculate the vector 'Hu_norm'
%
[Hu_norm] = setupHunorm(N, Hu_norm, th, h);

boundIndices = sort(boundaryUIndices(N));
allFluxes = 1:calculateAmountOfFluxes(N);

% Invert it. The boundaries should be thrown away
keepIndices = setdiff(allFluxes, boundIndices);

Hu_norm2 = H1t1(boundIndices, boundIndices);

%
%  Remove corresponding row and columns from the Hodge matrix and also
%  remove the corresponding 'rows' from Hu_norm
%
H1t1 = H1t1(keepIndices,keepIndices);
Hu_norm = Hu_norm(Hu_norm~=0); % Belonging to H1t1, and should be inversed for Ht11


%%
% Set up the incidence E^{2,1} between 1-cochain circulation and 2-cochain vorticity on
% the inner-oriented (extended) grid
%
% This incidence matrix will be called 'E21' in the program
%
[E21, totalExtendedFluxes] = setupE21(N);
fE21 = full(E21);

% Inserting prescribed tangential bundary conditions
% .................... Done, see below


%
% Remove columns from the incidence matrix E21 corresponding to both the
% prescribed tangental velocities and normal velocities
% .................... Done, see below


%
% Store the prescribed normal and tangential velocities in the vector
% 'u_pres'
[u_pres, E21] = extractUPresAndStripE21( E21, U_wall_top, U_wall_bot, U_wall_left, U_wall_right, ...
                                            V_wall_top, V_wall_bot, V_wall_left, V_wall_right, h );
fE21_final = full(E21);
fu_pres = full(u_pres);

%
% Set up the Hodge matrix which maps inner-oriented 2-cochains to
% outer-oriented 0-cochains. Assume that the vorticity is constant in the
% inner-oriented 2-cells. This will give a diagonal Hodge matrix. Call this
% Hodge matrix 'Ht02'
[Ht02] = setupHt02(N, h);



% Set up the Hodge matrix which maps inner-oriented 1-cochains to
% outer-oriented 1-cochains. Call this Hodge matrix 'Ht11'. Assume again
% that everything is constant over the 1-cells, which will then yield a
% diagonal Hodge matrix.
Ht11 = inv(H1t1);

    
%
% The prescribed velocties will play a role in the momentum equation
%
u_pres_orig = u_pres;
u_pres = H1t1*E21'*Ht02*u_pres;


%
% Now all matrices are set up and the time stepping can start. 'iter' will
% record the number of time steps. This allows you to give output after a
% preselected number of time steps.
%
% 'diff' will be the maximal du/dt or dv/dt. If 'diff' is sufficiently
% small, steady state has been reached. Determine a suitable value for
% 'tol'
%

diff = 1;
iter = 1;

%% Prepare convective term matrices
% Get incidence matrices for the convection
[ M1, M2, M3, M4, M5, M6, M7 ] = setupConvectionIncidenceMatrices( N );

% Get the indices of the prescribed fluxes
prescribedFluxIndices = extendedBoundaryUIndices(N);

% From the prescribed indices, get the indices of the unknown fluxes
unknownFluxesIndices = setdiff(1:totalExtendedFluxes, prescribedFluxIndices);

% Get the vector containing the prescribed fluxes (in order as described in
% this code)
prescribedFluxes = setupupres(N, U_wall_top, U_wall_bot, U_wall_left, U_wall_right, ...
                                            V_wall_top, V_wall_bot, V_wall_left, V_wall_right, h);
                                        
% Create the containing vector for all the fluxes. In this vector the
% unknown fluxes need to be updated every iteration.
uForIteration = ones(totalExtendedFluxes, 1);

% Fill in the prescribed fluxes
uForIteration(prescribedFluxIndices) = prescribedFluxes;

% Set up the matrix for the Poisson equation 
A = -tE21*(H1t1\tE21');       % from max


%% Start solver
Time.Matrices = toc;
disp(['Matrix creation time: ',num2str(Time.Matrices)])

iterationTimeKeeper=[];
iterationCounter=[];
iterationDiffKeeper=[];
while diff > tol
   
    % Fill in the unknown fluxes
    uForIteration(unknownFluxesIndices) = u;
    
    xi = E21 * u + u_pres_orig;
    txi = Ht02*xi;
    
    %
    %   Calculate the convective terms using the vorticty and the local
    %   velocity field. Store the convective terms in a vector
    %   'convective'.
    %
    %  Note that you only have to set up the convective terms for true
    %  velocity unknowns and not for those inner-oriented circulations for
    %  which the value is already known.
    %
    
    % - ([7]*{h})./(4*[1]*{h}) .* ([3]*{u}) .* ([4]*{xi}) -
    % ([7]*{h})./(4*[2]*{h}) .* ([5]*{u}) .* ([6]*{xi})
    
    convective =  - (M7*h')./(4*M1*h') .* (M3*uForIteration) .* (M4*txi)  ... % txi instead of xi from the assignment
                  - (M7*h')./(4*M2*h') .* (M5*uForIteration) .* (M6*txi);
    
    
    
    % Set up the right hand side for the Poisson equation for the pressure
    
    rhs_Poisson  =   tE21*(H1t1\(u/dt  - convective - H1t1*E21'*Ht02*E21*u/Re - u_pres/Re)) + u_norm/dt; 

    
    % Solve for the new pressure
    
    p = A\rhs_Poisson;
    
    % Store the velocity from the previous time step in the vector u_old
    
    uold = u;
    
    % Update the velocity field

      u = u - dt* (convective - tE21'*p + H1t1*E21'*Ht02*E21*u/Re + u_pres/Re); 

    %
    %  Every other 1000 iterations check whether you approach steady state
    %  and check whether you satisfy conservation of mass. The largest
    %  rate at which mass is destroyed or created is denoted by 'maxdiv'.
    %  This number should be very small, in the order of machine precision.
    
    if mod(iter,1000) == 0
    
        maxdiv = max(tE21*(H1t1\u) + u_norm);
        
        diff = max(abs(u-uold))/dt;
        
        % Log
        iterationCounter = [iterationCounter;iter];
        iterationDiffKeeper = [iterationDiffKeeper;diff];
        iterationTimeKeeper = [iterationTimeKeeper;toc];
        elapsedTime = iterationTimeKeeper(end);
        
        disp(['Iteration: ',num2str(iter),'    maxDiv: ',num2str(maxdiv),'   diff: ',num2str(diff)]);
    end
        
       
    iter = iter + 1;
    
end
Time.Iteration = toc - Time.Matrices;
disp(['Iteration time: ',num2str(Time.Iteration)]);

%% 
%   Save all variables to .mat file 
%
saveFolder = './Results/';
fileName = sprintf('%s%g%s%g%s%s%s%s%s','N',N,'Tol',tol,'_',datestr(now,'yyyymmdd_HHMM'),'_',getenv('USERNAME'),'.mat');
savepath = [saveFolder,fileName];
save(savepath);
%%
%  Produce desired out put to compare your results with those given in the
%  reference by Botella and Peyret
%

res = Result.load(savepath);

res.plotConvergence(1,1)
res.plotConvergence(2,1)
res.plotConvergence(3,1)

res.plotStreamLines()

res.plotVorticity(2,1)

res.plotPressure()

res.plotValidation(1,'h',1)
res.plotValidation(2,'h',1)
res.plotValidation(3,'h',1)
res.plotValidation(1,'v',1)
res.plotValidation(2,'v',1)
res.plotValidation(3,'v',1)

Time.PostProcessing = toc - Time.Iteration;

disp(['PostProcessing time: ',num2str(Time.PostProcessing)]);