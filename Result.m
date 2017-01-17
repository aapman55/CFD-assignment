%% Assignment for AE4134: CFD I
% Created by:
% Zhi-li Liu 4146557
% Jasper van Wensveen 4142179
%
%RESULT This class contains the result matrix of the simulation. This
%class is made for ease of plotting and reducing code duplication.
classdef Result < handle

    
    properties
        res;                    % The result matrix
        fontsize;               % Standard fontsize for all plots
        stdFigSize;             % Standard Figure sizes for all plots
        avgU;                   % The average u velocity at the inner 0 cells
        avgV;                   % The average v velocity at the inner 0 cells
        statP;                  % The static pressure at the inner 0 cells
        matAvgU;                % The average u velocity vector turned into matrix format
        matAvgV;                % The average v velocity vector turned into matrix format
        vortIsoLines;           % The ISO lines used in Botella and Peyret
        presIsoLines;           % The ISO lines used in Botella and Peyret
        streamIsoLines;         % The ISO lines used in Botella and Peyret
    end
    
    methods
        function obj = Result(resultMatrix, fontsize, stdFiguresize)
           obj.res = resultMatrix; 
           obj.fontsize = fontsize;
           obj.stdFigSize = stdFiguresize;
           
           % Pre processing
           r = obj.res;
           % Get the incidence matrix for the average of U and V at the pressure
            % points
            [incAvgU, incAvgV] = setupIncidenceMatrixAvgUV(r.N);

            % Create full inner circulation vector including prescribed circulation
            % (because the calculated u's are on the inner grid)
            [fullInnerOrientedCirc, unknownUIndices] = createFullInnerCirculationVector(r.N, r.U_wall_left, r.U_wall_right, r.V_wall_bot, r.V_wall_top);
            fullInnerOrientedCirc(unknownUIndices) =  r.u;

            % Set-up the array containing the circulation (inner-grid) lengths for each velocity
            outerHeights = [repmat(r.h,1,r.N),repelem(r.h,r.N)];

            % Turn the circulations into velocities (velocities) by dividing the circulations by
            % the length of the segment through which the circuulation goes

            fullInnerOrientedCirc = fullInnerOrientedCirc./outerHeights';

            uCirc = fullInnerOrientedCirc(1:end/2);
            vCirc = fullInnerOrientedCirc(end/2+1:end);

            % Get the actual values for the average u and v
            obj.avgU = incAvgU*uCirc;
            obj.avgV = incAvgV*vCirc;            
                        
            % Put avgU and avgV in matrix format
            obj.matAvgU = assignmentReshape(obj.avgU,r.N,r.N);
            obj.matAvgV = assignmentReshape(obj.avgV,r.N,r.N);
            
            % statP no Hodge
            obj.statP = r.p - 0.5.*(obj.avgU.^2 + obj.avgV.^2); 
            
            % Stream function ISO lines
            obj.streamIsoLines = [  0.1175, 0.115, 0.11, 0.1, 9E-2, 7E-2, 5E-2, 3E-2, 1E-2,...
                                    0, -1E-6, -1E-5, -5E-5, -1E-4, -2.5E-4, -5E-4, -1E-3, -1.5E-3];
        end
%% ====================
% Pressure plot
%======================     
        function plotPressure(obj)
            % Shorten result notation to r
            r = obj.res;
                        
            % Create hodge matrix H0t2
            % H0t2 = setupH0t2(th);

            % Calculate the real static pressure (p). P = p + 0.5*(u^2 + v^2) 
            % statPHodge = H0t2*p - 0.5.*(avgU.^2 + avgV.^2); 

            % The x locations are the locations of the inner grid without the extended
            % boundary conditions (So removing the first and last entry)
            X = r.x(2:end-1);

            % Create the 2D grid
            [X1,Y1] = meshgrid(X,X);

            figure('position',[0 0 obj.stdFigSize])
            % pcolor(X1,Y1,reshape(p-round(mean(p)),N,N))
            % pcolor(X1,Y1,assignmentReshape(statP,N,N))
            % contourf(X1,Y1,assignmentReshape(statPHodge,N,N),linspace(-1.5E3,1E3,50))

            % Try with no Hodge p only
            contourf(X1,Y1,assignmentReshape((obj.statP),r.N,r.N),linspace(-.8,.2,100))
            axis image
            ylabel(colorbar,'Pressure [N/m^2]')
            xlabel('Horizontal direction [m]')
            ylabel('Vertical direction [m]')
            set(gca,'fontsize',obj.fontsize) 
        end
        
%% ===================================
% Stream lines (poging tot)
%=====================================
% The streamfunction at a certain point is defined as the integral of (u dy
% - v dx) integrated from the datum point to that certain point (ds). As the path
% taken does not matter dx will be considered as the x difference between that poit and the datum 
%, dy the difference in y location of that point with the datum and ds the length calculated with
% Pythagoras with dx and dy.

        function plotStreamLines(obj)
            r = obj.res;
            
            % The x locations are the locations of the inner grid without the extended
            % boundary conditions (So removing the first and last entry)
            X = r.x(2:end-1);

            % Create the 2D grid
            [X1,Y1] = meshgrid(X,X);
            % The datum point is taken to be [0.5, 0.5]
            dX1ForStreamline = X1-0.47;
            dY1ForStreamline = Y1-0.6;

            streamFunction = (obj.matAvgU.*dY1ForStreamline - obj.matAvgV.*dX1ForStreamline).*sqrt(dX1ForStreamline.^2+dY1ForStreamline.^2);

            figure('position',[0 0 obj.stdFigSize])
            contourf(X1,Y1,streamFunction,linspace(-0.2,0.1,75))
            axis image
            ylabel(colorbar,'Streamfunction')
            xlabel('Horizontal direction [m]')
            ylabel('Vertical direction [m]')
            set(gca,'fontsize',obj.fontsize)
        end
%% ==============================
% Plot the vorticity XI
%================================        
        function plotVorticity(obj, variant, createFig)
            r = obj.res;
            
            % Turn the 1D vector txi into the 2D matrix corresponding to the numbering
            % scheme as used in the assignment
            XI = assignmentReshape(r.txi,r.N+1,r.N+1);
            X = r.x(1:end-1)+r.h/2;
            [X2,Y2]=meshgrid(X,X);
            
            %Create new figure
            if exist('createFig','var') && createFig
                figure('position',[0 0 obj.stdFigSize])
            end
            
            switch variant
                case 1         
                    pcolor(X2,Y2,XI)
                case 2
                    contourf(X2,Y2,XI,-30:15)
            end
            
            % If no new figure is created, do not change figure settings
            if exist('createFig','var') && createFig
                caxis([-30,15])
                axis image
                ylabel(colorbar,'Vorticity [1/s]')
                xlabel('Horizontal direction [m]')
                ylabel('Vertical direction [m]')
                set(gca,'fontsize',obj.fontsize)
            end
            
        end
%% =========================
%   Make convergence relevant plots
% ==========================        
        function plotConvergence(obj, variant, createFig)
            
            r = obj.res;
            if exist('createFig','var') && createFig
                figure('position',[0 0 obj.stdFigSize])
            end

            switch variant
                case 1                    
                    plot(r.iterationCounter, r.iterationDiffKeeper);
                    xlabel('Iteration [-]')
                    ylabel('Maximum change in u [m^2/s]')
                case 2
                    loglog(r.iterationCounter, r.iterationDiffKeeper);
                    xlabel('Iteration [-]')
                    ylabel('Maximum change in u [m^2/s]')
                case 3
                    loglog(r.iterationTimeKeeper, r.iterationDiffKeeper);
                    xlabel('Time [s]')
                    ylabel('Maximum change in u [m^2/s]')
                case 4                    
                    plot(r.iterationTimeKeeper, r.iterationDiffKeeper);
                    xlabel('Time [s]')
                    ylabel('Maximum change in u [m^2/s]')
                case 5                    
                    plot(r.iterationCounter, r.iterationTimeKeeper);
                    xlabel('Iterations [-]')
                    ylabel('Time [s]')
                case 6                    
                    loglog(r.iterationCounter, r.iterationTimeKeeper);
                    xlabel('Iterations [-]')
                    ylabel('Time [s]')
                otherwise
                    warning('This plot variant does not exist!')
            end
            
            if exist('createFig','var') && createFig
                grid minor
                axis tight
                set(gca,'fontsize',obj.fontsize)
            end
        end
%% =========================
% Validation
% ==========================     
% Variant 1 = pressure, 2 = vorticity, 3 = velocity
        function plotValidation(obj, variant, HV, createFig)
            r = obj. res;
            
            % Input checking
            if(variant ~= 1 && variant ~= 2 && variant~=3)
                error('The variant can only be 1,2 or 3!')
            end
            
            if(HV ~= 'h' && HV~= 'v')
                error('The HV can only be h or v')
            end
            
            % Load validationData
            horizontalCenterlineData = dlmread('validationHorCenterline',';',1,0);
            verticalCenterlineData = dlmread('validationVerCenterline',';',1,0);
            
            % Create new variable for make Fig
            makeFig = false;
            
            % Create optional figure
            if exist('createFig','var') && createFig
                figure('position',[0 0 obj.stdFigSize])
                makeFig = true;
            end
            
            % vertical velocity
            Hv = obj.matAvgV(end/2,:);
            
            % horizontal velocity
            Vu = obj.matAvgU(:,end/2);
                        
            switch variant
                case 1
                %============================
                % Pressure
                %============================

                % Put pressure in matrix format
                matStatPres = assignmentReshape(obj.statP,r.N,r.N);
                
                switch HV
                    case 'h'
                        % Pressure
                        Hpres = matStatPres(end/2,:);
                        % Own simulation
                        if(makeFig)
                            % Botella validation data
                            plot(horizontalCenterlineData(:,1),horizontalCenterlineData(:,4),'k--','linewidth',2)
                        end
                        hold on
                        plot(r.x(2:end-1),Hpres)

                        % Pressure shifted with a constant to match botella data
                        plot(r.x(2:end-1),Hpres+horizontalCenterlineData(1,4)-Hpres(1))
                        xlabel('Horizontal direction [m]')
                    case 'v'
                        % Pressure
                        Vpres = matStatPres(:,end/2);
                        % Own simulation
                        if(makeFig)
                            % Botella validation data
                            plot(verticalCenterlineData(:,1),verticalCenterlineData(:,4),'k--','linewidth',2)
                        end
                        hold on
                        plot(r.x(2:end-1),Vpres)

                        % Pressure shifted with a constant to match botella data
                        plot(r.x(2:end-1),Vpres+verticalCenterlineData(end,4)-Vpres(1))
                        xlabel('Vertical direction [m]')
                end
                if(makeFig)
                    ylabel('Pressure [N/m^2]')
                    legend('Botella and Peyret Validation','Own simulation','Own simulation with constant','location','North')
                    grid minor
                    set(gca,'fontsize',obj.fontsize)
                end
                
                case 2
                    %============================
                    % Vorticity
                    %============================
                    % Turn the 1D vector txi into the 2D matrix corresponding to the numbering
                    % scheme as used in the assignment
                    XI = assignmentReshape(r.txi,r.N+1,r.N+1);
                    X = r.x(1:end-1)+r.h/2;
                    [X2,Y2]=meshgrid(X,X);
                    
                    switch HV
                        case 'h'
                            % Plot own simulation vorticity data
                            if(makeFig)
                                % Plot botella vorticity data
                                plot(horizontalCenterlineData(:,1),horizontalCenterlineData(:,5),'k--','linewidth',2)
                            end
                            hold on
                            plot(X2(floor(end/2),:), XI(floor(end/2),:))
                            xlabel('Horizontal direction [m]')

                        case 'v'
                            % Plot own simulation vorticity data
                            if(makeFig)
                                % Plot botella vorticity data
                                plot(verticalCenterlineData(:,1),verticalCenterlineData(:,5),'k--','linewidth',2)
                            end
                            hold on
                            plot(Y2(:,floor(end/2)), XI(:,floor(end/2)))
                            xlabel('Vertical direction [m]')
                    end
                    if(makeFig)
                        ylabel('Vorticity [1/s]')
                        legend('Botella and Peyret Validation','Own simulation','location','south')
                        grid minor
                        set(gca,'fontsize',obj.fontsize)
                    end
                    
                case 3
                    %============================
                    % Velocity
                    %============================
                    % The x locations are the locations of the inner grid without the extended
                    % boundary conditions (So removing the first and last entry)
                    X = r.x(2:end-1);

                    % Create the 2D grid
                    [X1,Y1] = meshgrid(X,X);
                    
                    switch HV
                        case 'h'
                            % Plot own simulation vertical velocities
                            if(makeFig)
                                plot(horizontalCenterlineData(:,1),horizontalCenterlineData(:,3),'k--','linewidth',2)
                            end
                            hold on;
                            plot(X1(end/2,:),Hv)
                            xlabel('Horizontal direction [m]')
                            ylabel('Vertical velocity [m/s]')
                        case 'v'                            
                            % Plot own simulation horizontal velocities
                            if(makeFig)
                                plot(verticalCenterlineData(:,1),verticalCenterlineData(:,3),'k--','linewidth',2)
                            end
                            hold on
                            plot(Y1(:,end/2),Vu)
                            xlabel('Vertical direction [m]')
                            ylabel('Horizontal velocity [m/s]')
                    end
                    legend('Botella and Peyret Validation','Own simulation','location','south')
                    grid minor
                    set(gca,'fontsize',obj.fontsize)
            end

        end
    end
    
    methods(Static)
        function resultObject = load(path)
            resultObject = Result(load(path), 16, [800, 600]);
        end
    end
    
end

