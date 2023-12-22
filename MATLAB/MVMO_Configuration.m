%**************************************************************************
%           MVMO-based classifier identification parameters
%**************************************************************************
%By: Jaime Cepeda
%Modified by: José Enríquez
% =========================================================================
%                         About MVMO Algorithm 
% =========================================================================
% Programmer: W. Nakawiro
% Constraint handling: Modified self-adaptive penalty scheme 
% Date:  3 Nov. 10
% Reference: 
% I. Erlich, G.K. Venayagamoorthy and W.Nakawiro,"A mean-variance
% optimization algorithm", In.Proc. 2010 IEEE World Congress on
% Computational Intelligence, Barcelona, Spain, July 2010
% =========================================================================


% =========================================================================
%                         Setting global variables 
% =========================================================================
% This global structure stores infomation avaliable to all functions
global parameter
global printff sn
global iopt % for predictive control optimization
% =========================================================================

% =========================================================================
%                     Objective Function Selection

% Select a number depending on the objective function you choose

% (1): tracking logic
% (2): rate limiting logic
% (3): ACE filtering logic

iss=3;

% % =========================================================================
% %                     Setting limits for optimization 
% % =========================================================================
% 
% % Lower Limits: LB_mtrx_est
% % Upper Limits: LU_mtrx_est
% 

switch iss
    case 1
        %[T1;T2]
        LB_mtrx_est=[4; 40];
        LU_mtrx_est=[500; 80]; 
        
    case 2
        %[T3;T4]
        LB_mtrx_est=[4; 10];
        LU_mtrx_est=[500; 40];
   
    case 3
        %[window size;XWMATC]
        LB_mtrx_est=[1;    4];
        LU_mtrx_est=[45;    100];
end

% =========================================================================


% =========================================================================
%                             Define MVO parameters
% =========================================================================

max_run = 1; % maximum number of runs 
printff= 10; % print  the results at every printff FEs on the screen
parameter.MaxEval =2000; % max no. of function/fitness evaluations
parameter.n_random_max = 2; %number of variables to vary
parameter.n_random_min = 4;
parameter.n_random=2; %N° variables 

%--------------------------------------------------------------------------
minTime = Inf;
parameter.mode =4; %variable selection method
parameter.n_tosave = 4;
parameter.n_random_freq = 1000; 
parameter.n_random_begin = 1500;

parameter.fs = 1.05;        %shaping scaling factor    
parameter.fsmin = 1.05;
parameter.fsmax = 1.05;

parameter.AF = 2.0d0;       %Assymetry factor
parameter.sd = 0.01;        %initial value of shape factor

sn = true; % logic control for screen printing
    
saveff =100;
nnsave =floor(parameter.MaxEval/saveff)+1; 
FConv = NaN*ones(nnsave,max_run);
%--------------------------------------------------------------------------
%Datos
parameter.n_var = length(LU_mtrx_est);   % total no. of variables 
parameter.ncon_var = length(LU_mtrx_est); %length(xx_yy);  % no. of continuous variables
parameter.ndis_var = parameter.n_var-parameter.ncon_var;   % no. of discrete variables
parameter.n_constr = 0; % no. of constraints

%--------------------------------------------------------------------------
       
% Min and max boundaries 
parameter. x_min = LB_mtrx_est'; %optimization variables (x,y)
parameter. x_max = LU_mtrx_est';

parameter. v_min = LB_mtrx_est'; %variance
parameter. v_max = LU_mtrx_est';
% =========================================================================


%==========================================================================    
%                 Performing MVMO-based classifier optimization
%========================================================================== 
XBest=zeros(max_run,parameter.n_var);
Iter_nume_best=zeros(parameter.MaxEval,max_run);
Objv_func_best=zeros(parameter.MaxEval,max_run);
Func_feas_best=zeros(parameter.MaxEval,max_run);

for iopt = 1:max_run
    Best_solu_best(iopt).param=zeros(parameter.MaxEval,parameter.n_var);
end
FConv_best=zeros(floor(parameter.MaxEval/saveff)+1,max_run);
Ctime_elap = zeros(1,max_run);
Ctime_min = zeros(1,max_run);

tic;
for iopt = 1:max_run
    tStart = tic;
    
    %Normalization   
    Scaling = (parameter.x_max-parameter.x_min);
    %you can generate N random numbers in the interval (a,b) with the formula r = a + (b-a).*rand(N,1)
    x_normalized = parameter.v_min+(parameter.v_max-parameter.v_min).*rand(1,parameter.n_var);
    %Xnorm=X-Xmin/Xmax-Xmin
    x_normalized = (x_normalized-parameter.x_min)./Scaling;

    %-----------------------------------------------------------------------
    % Call objective function
    [FConv,xbest_norm,Iter_nume, Objv_func, Func_feas, Best_solu] = MVMO_old_nakawiro(x_normalized,iopt,saveff,iss,parameter);
    %-----------------------------------------------------------------------
    
    tElapsed = toc(tStart);
    Ctime_elap(iopt)=tElapsed;
    Ctime_min(iopt)=min(tElapsed, minTime);
        
    XBest(iopt,:) = parameter.x_min+(parameter.x_max-parameter.x_min).*xbest_norm;
    Iter_nume_best(:,iopt)=Iter_nume; 
    Objv_func_best(:,iopt)=Objv_func; 
    Func_feas_best(:,iopt)=Func_feas;
    Best_solu_best(iopt).param=Best_solu;
    FConv_best(:,iopt)=FConv;
end
averageTime = toc/iopt;
%==========================================================================   


%==========================================================================    
%                       Save results and plot
%==========================================================================
results.Ctime_elap=Ctime_elap;
results.Ctime_min=Ctime_min;
results.Iter_nume_best=Iter_nume_best;
results.Objv_func_best=Objv_func_best;
results.Func_feas_best=Func_feas_best;
results.Best_solu_best=Best_solu_best;
results.XBest=XBest;
results.FConv_best=FConv_best;

save results_mvmo_erlich  results averageTime iopt 

%figure iter number vs objective function evaluation
figure
plot(results.Iter_nume_best,results.Objv_func_best)
xlabel('Iteration number')
ylabel('Objective function evaluation')

%figure iter number vs x1 
figure
plot(results.Iter_nume_best,results.Best_solu_best.param(:,1))
xlabel('Iteration number')
ylabel('x1')
%figure iter number vs x2 
figure
plot(results.Iter_nume_best,results.Best_solu_best.param(:,2))
xlabel('Iteration number')
ylabel('x2')

% save Opt_mean_param_mvmo Opt_mean_param
fprintf('\nx1_opt = %.4f\n',results.XBest(1))
fprintf('x2_opt = %.4f\n',results.XBest(2))






