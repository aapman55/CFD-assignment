%% Assignment for AE4134: CFD I
% Created by:
% Zhi-li Liu 4146557
% Jasper van Wensveen 4142179
%% Results comparison (between N, tolerances etc.
clear; close all; clc;

stdFigureSize = [800, 600];     % Width and height
fontsize = 16;
%% ==========================
% Tolerances (N = 32)
%============================
N32Tol(1) = Result.load('./Results/correcteResults/N32Tol0.01_20161123_0932_zl.mat');
N32Tol(2) = Result.load('./Results/correcteResults/N32Tol0.001_20161116_2233_Jasper.mat');
N32Tol(3) = Result.load('./Results/correcteResults/N32Tol0.0001_20161116_2235_Jasper.mat');
N32Tol(4) = Result.load('./Results/correcteResults/N32Tol1e-05_20161116_2238_Jasper.mat');
N32Tol(5) = Result.load('./Results/correcteResults/N32Tol1e-06_20161116_2242_Jasper.mat');
N32Tol(6) = Result.load('./Results/correcteResults/N32Tol1e-07_20161116_2248_Jasper.mat');


for i = N32Tol   
    i.plotVorticity(2,1);    
end

N32Tol(1).plotValidation(3,'v',1)
N32Tol(2).plotValidation(3,'v',0)
N32Tol(3).plotValidation(3,'v',0)
N32Tol(4).plotValidation(3,'v',0)
N32Tol(5).plotValidation(3,'v',0)
N32Tol(6).plotValidation(3,'v',0)

legend( 'Botella and Peyret',...
        ['N=32, Tol = ',num2str(N32Tol(1).res.tol)],...
        ['N=32, Tol = ',num2str(N32Tol(2).res.tol)],...
        ['N=32, Tol = ',num2str(N32Tol(3).res.tol)],...
        ['N=32, Tol = ',num2str(N32Tol(4).res.tol)],...
        ['N=32, Tol = ',num2str(N32Tol(5).res.tol)],...
        ['N=32, Tol = ',num2str(N32Tol(6).res.tol)])
grid minor

%% ==========================
% Turbo Boost (N=32 , Tol 1E-5)
%============================
Turbo(1) = Result.load('./Results/correcteResults/N32Tol1e-05_20161116_2214_zhililiu.mat');
Turbo(2) = Result.load('./Results/correcteResults/N32Tol1e-05_20161116_2203_zhililiu.mat');
Turbo(3) = Result.load('./Results/correcteResults/N32Tol1e-05_20161116_2206_zhililiu.mat');
Turbo(4) = Result.load('./Results/correcteResults/N32Tol1e-05_20161116_2207_zhililiu.mat');
Turbo(5) = Result.load('./Results/correcteResults/N32Tol1e-05_20161116_2208_zhililiu.mat');

% Make plots
Turbo(1).plotConvergence(3,1); hold on;
Turbo(2).plotConvergence(3,0)
Turbo(3).plotConvergence(3,0)
Turbo(4).plotConvergence(3,0)
Turbo(5).plotConvergence(3,0)
legend(['dt = ',num2str(Turbo(1).res.dt)],...
        ['dt = ',num2str(Turbo(2).res.dt)],...
        ['dt = ',num2str(Turbo(3).res.dt)],...
        ['dt = ',num2str(Turbo(4).res.dt)],...
        ['dt = ',num2str(Turbo(5).res.dt)])
    
% Make more plots    
Turbo(1).plotConvergence(4,1); hold on;
Turbo(2).plotConvergence(4,0)
Turbo(3).plotConvergence(4,0)
Turbo(4).plotConvergence(4,0)
Turbo(5).plotConvergence(4,0)
legend(['dt = ',num2str(Turbo(1).res.dt)],...
        ['dt = ',num2str(Turbo(2).res.dt)],...
        ['dt = ',num2str(Turbo(3).res.dt)],...
        ['dt = ',num2str(Turbo(4).res.dt)],...
        ['dt = ',num2str(Turbo(5).res.dt)])
    
% Trend plot
times = [Turbo(1).res.dt,Turbo(2).res.dt,Turbo(3).res.dt,Turbo(4).res.dt,Turbo(5).res.dt];
timekeeper = [  Turbo(1).res.iterationTimeKeeper(end),...
                Turbo(2).res.iterationTimeKeeper(end),...
                Turbo(3).res.iterationTimeKeeper(end),...
                Turbo(4).res.iterationTimeKeeper(end),...
                Turbo(5).res.iterationTimeKeeper(end)];
loglog([1,2,3,4,5],timekeeper)
grid minor
axis tight
%% ==========================
% Effect of N Tol 1E-6
%============================    
Ncompare(1) = Result.load('Results/correcteResults/N8Tol1e-06_20161116_2231_Jasper.mat');
Ncompare(2) = Result.load('Results/correcteResults/N16Tol1e-06_20161116_2232_Jasper.mat');
Ncompare(3) = Result.load('Results/correcteResults/N32Tol1e-06_20161116_2242_Jasper.mat');
Ncompare(4) = Result.load('Results/correcteResults/N48Tol1e-06_20161122_2033_Jasper.mat');
Ncompare(5) = Result.load('Results/correcteResults/N56Tol1e-06_20161124_1013_Jasper.mat');
Ncompare(6) = Result.load('Results/correcteResults/N64Tol1e-06_20161117_1139_Jasper.mat');

for i = Ncompare
   i.plotVorticity(2,1)
end

% Convergence
Ncompare(1).plotConvergence(3,1);hold on;
Ncompare(2).plotConvergence(3,0);
Ncompare(3).plotConvergence(3,0);
Ncompare(4).plotConvergence(3,0);
Ncompare(5).plotConvergence(3,0);
Ncompare(6).plotConvergence(3,0);
legend(['N = ',num2str(Ncompare(1).res.N)],...
        ['N = ',num2str(Ncompare(2).res.N)],...
        ['N = ',num2str(Ncompare(3).res.N)],...
        ['N = ',num2str(Ncompare(4).res.N)],...
        ['N = ',num2str(Ncompare(5).res.N)],...
        ['N = ',num2str(Ncompare(6).res.N)])
    
Ncompare(1).plotConvergence(6,1);hold on;Ncompare(2).plotConvergence(6,0);
Ncompare(3).plotConvergence(6,0);
Ncompare(4).plotConvergence(6,0);
Ncompare(5).plotConvergence(6,0);
Ncompare(6).plotConvergence(6,0);
legend(['N = ',num2str(Ncompare(1).res.N)],...
        ['N = ',num2str(Ncompare(2).res.N)],...
        ['N = ',num2str(Ncompare(3).res.N)],...
        ['N = ',num2str(Ncompare(4).res.N)],...
        ['N = ',num2str(Ncompare(5).res.N)],...
        ['N = ',num2str(Ncompare(6).res.N)])

 % Needed time trend   

 x = zeros(6,1);
 y = zeros(6,1);
 for i = 1:length(Ncompare)
     x(i) = Ncompare(i).res.N;
     y(i) = Ncompare(i).res.iterationTimeKeeper(end);
 end
 
 fit(x,y,'exp1');
 
 %Manually read off the fit
 extraFit = @(x) 19.77*exp(0.09937*x);
 
 figure('position',[0 0 Ncompare(1).stdFigSize])
 loglog(x,y,'-s')
 hold on
 loglog(128,extraFit(128),'s')
 grid minor
 xlabel('N [-]')
 ylabel('Time [s]')
 set(gca,'fontsize',Ncompare(1).fontsize)
 

%% ===========================
% Validation plots
%=============================
N64Tol6 = Result.load('Results/correcteResults/N64Tol1e-06_20161117_1139_Jasper.mat');

N64Tol6.plotStreamLines()
N64Tol6.plotVorticity(2,1)
N64Tol6.plotPressure()

N64Tol6.plotValidation(1,'h',1)
N64Tol6.plotValidation(2,'h',1)
N64Tol6.plotValidation(3,'h',1)
N64Tol6.plotValidation(1,'v',1)
N64Tol6.plotValidation(2,'v',1)
N64Tol6.plotValidation(3,'v',1)