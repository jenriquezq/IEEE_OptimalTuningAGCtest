function [Fx1,xbest,Iter_nume, Objv_func, Func_feas, Best_solu] = MVMO_old_nakawiro(x_normalized,krun,saveff,iss,varargin)

global parameter;
global x_normalized_best  weight
global considered changed fs_factor
global printff sn;

iter_count1=0;

% Intialize the control parameters
n_var = parameter.n_var; % Number of problem variables
max_eval = parameter.MaxEval; % Max. number of fitness evaluations

% ----- data structure for the table ------
table.bests = zeros(parameter.n_tosave,parameter.n_var);
table.objective = zeros(parameter.n_tosave,1);
table.fitness = 1e10*ones(parameter.n_tosave,1);
table.feasibility = zeros(parameter.n_tosave,1);
%------------------------------------------

Iter_nume=zeros(max_eval,1);
Objv_func=zeros(max_eval,1);
Func_feas=zeros(max_eval,1);
Best_solu=zeros(max_eval,n_var);
Best_solu_norm=zeros(max_eval,n_var);


nsave = round(parameter.MaxEval/printff)+1;
insave = 0;
mean_save = zeros(nsave,parameter.n_random);
idx_save=mean_save; 
shape_save1 = zeros(nsave,parameter.n_random);
shape_save2=shape_save1;

mode =parameter.mode;
dddd=parameter.sd; %initial value of shape factor
delta_ddd= 0.0505d0/real(n_var)+1.d0; %Initial value of alternative shape
asymmetry=parameter.AF;

Fx_best =1.0d50; %1e50

n_randomly=parameter.n_random; % no. of variables that are changed randomly
n_to_save =parameter.n_tosave; % no. of solutions to be saved in the table
fs_factor = parameter.fs;

% Declaring arrays
considered = true(1,n_var);
variance =ones(1,n_var);
vsel = 1:n_randomly;        %vsel vector Number of variable changed randomly
vv =1:n_randomly;       %vv vector with Number of variable changed randomly
no_in = 0;      %Initial number of variables selected for mutation
no_inin =0;

if n_randomly > n_var
    n_randomly=n_var;
end

x_normalized_best = x_normalized;
meann = x_normalized;

kidx =0;
for k=1:max_eval
    nnn = 1000; nnn1 = 6000; nnp=2;
    iter_count1=iter_count1+1;
    
    if k<=nnn
        n_randomly =parameter.n_random;
    elseif (k>nnn)&& (k<=nnn1)
        mslope = (parameter.n_random-nnp)/(nnn-nnn1);
        n_randomly =  round(mslope*(k-nnn)+parameter.n_random);%round(parameter.n_random-k*(parameter.n_random-ntapp)/(nnn1-nnn));
    elseif (k> nnn1)
        n_randomly =nnp;
    end
    
     if n_randomly <nnp;  n_randomly = nnp; end
    % Evaluate fitness of the individual

%----------------------------------------------------------------------------------------------------
%OBJECTIVE FUNCTION
    %ffx function objective
    %ggx inequality constraint
      [ffx, ggx,x_normalized] = Objective_Function(x_normalized,iss);
%----------------------------------------------------------------------------------------------------    
           
        xt.objective = ffx; %set result for objective function
        xt.fitness = static_penalty(ffx,ggx); 
        if xt.fitness==xt.objective
            xt.feasibility= true;
        else
            xt.feasibility= false;
        end
    % apply the rule-based elitism 
    isxbetter = elitism(xt,table,1);
    if isxbetter
         Fx_best = xt.objective; %Make Fx_best equal to value of objective function
    end
    % storing the n-best solutions to the archive
    save_best_feasible();
    % printing the global best result to the screen
    if ((k==1)||(mod(k,printff) ==0))&& sn       
            if sn
             fprintf('%7d %7d %.2e %.2e  %7d  \n',...
             krun, k,table.objective(1),table.fitness(1),table.feasibility(1));               
            end
    end

    if ((k==1)||(mod(k,saveff) ==0))
        kidx = kidx+1;
        Fx1(kidx,:) = Fx_best;
    end
    
        x_normalized = x_normalized_best;

        considered(1:n_var) = false;
        % call the random variable selection strategies
        VariableSelect1();%,n_var,n_randomly);
        % change the corresponding variable(s)
        vrand = rand(1,n_var);
        x_normalized(considered) = vrand(considered);
    isvv =0;
    for ivar=1:n_var
        if considered(ivar)
            sss1 = -log(variance(ivar))*fs_factor ; %shape variable calculation, variance is initialized by one
            sss2=sss1;
            
            if x_normalized_best(ivar) < meann(ivar)  %x_normalized_best is equal meann for the first time 
                   sss2 = asymmetry*sss2 ;
            elseif x_normalized_best(ivar) > meann(ivar) 
                   sss1=asymmetry*sss1;
            else
                if sss2 >= dddd
                    dddd = dddd*delta_ddd;
                else
                    dddd = dddd/delta_ddd;
                end
                sss1 = dddd;
            end
             x_normalized(ivar)=  h_function(meann(ivar),sss1,sss2,x_normalized(ivar));
              
             if (~mod(k,printff))||(k==1)
                 if isvv ==0
                 insave = insave+1;
                 end
                 isvv = isvv+1;
                 mean_save(insave,isvv) = meann(ivar);
                 shape_save1(insave,isvv) = sss1;
                 shape_save2(insave,isvv) = sss2;
                 idx_save(insave,isvv) = ivar;
             end
             
        end
    end  
    
    Iter_nume(iter_count1,:)=k;
    Objv_func(iter_count1,:)=table.objective(1,1);
    Func_feas(iter_count1,:)=table.feasibility(1,1);
    Best_solu(iter_count1,:)=parameter.x_min+(parameter.x_max-parameter.x_min).*table.bests(1,:); %denormalize
    Best_solu_norm(iter_count1,:)=table.bests(1,:);
end
xbest = table.bests(1,:);
%%
    function save_best_feasible()
    % store best
    persistent ixbetter;
    no_in = no_in+1;
    changed = false;
    ixbetter = false;
   
    if no_in ==1 % the solution coming to the table for the first time
        table.bests(1,:) = x_normalized;
        table.objective(1,:) = xt.objective;
        table.fitness(1,:) = xt.fitness;
        table.feasibility(1,:) =xt.feasibility;
        
         x_normalized_best = table.bests(1,:); % set the best solution to the one of the first rank
         weight =1.d0/double(n_to_save); % weight for mean calculation
         no_inin=no_inin+1;
    else % not for the first time and check for the update
       i_position =0;
       % check if the new coming solution is less than any in the table
        for i=1:n_to_save 
            if changed; break; end
            ixbetter = elitism(xt,table,i);
            if ixbetter %Fx <= table.objective(i,1) 
                i_position = i; %if the new solution is better than any of the previous ones, it ocuppies the position of that solution
                changed =true;
                no_inin = no_inin+1; % how many times good solutions are in    
            end
        end
    end
    
    if changed % if the new individual is better than any archived individual
            % Move the individuals and corresponding fitness values
            % downward so that the individuals are sorted based on the
            % fitness value in a descending order 
            if (n_to_save >= i_position+1)
                for i=n_to_save:-1:i_position+1
                    table.objective(i,:) = table.objective(i-1,:);
                    table.bests(i,:) = table.bests(i-1,:);
                    table.fitness(i,:)=table.fitness(i-1,:);
                    table.feasibility(i,:)= table.feasibility(i-1,:);
                end
            end
            % save the new best
            table.bests(i_position,:) = x_normalized;
            table.objective(i_position,:) = xt.objective;
            table.fitness(i_position,:) = xt.fitness;
            table.feasibility(i_position,:) =xt.feasibility;
         
            % Set new global best
            x_normalized_best = table.bests(1,:);
            
         % calculation of mean and variance
         if (no_inin >= n_to_save)
             
                 if(i_position ~= 1)
                     idcontrol = true(1,n_var)|considered;
                 else
                     idcontrol = false(1,n_var)|considered;
                 end
                 mean_temp = weight*sum(table.bests(1:n_to_save,:),1);
                 meann( idcontrol) = mean_temp( idcontrol);
                 dummy_ctrl=100*ones(length(mean_temp),1);
                 dummy_ctrl( idcontrol) = mean_temp( idcontrol);
                 contro_diff(1,:)=meann-mean_temp;
                 contro_diff(2,:)=idcontrol;
                 contro_diff(3,:)=dummy_ctrl;
                 variance_temp = weight*sum((table.bests(1:n_to_save,:)-repmat(mean_temp,n_to_save,1)).^2,1);
                 id_nonzero = (variance_temp > 1.d-200)&( idcontrol);
                 variance(id_nonzero) = variance_temp( id_nonzero);  
         end

         
    end  % end "changed" 
    end % end "save_best_feasible"


    function VariableSelect1()
       if ismember(mode,[2,3,4]);
           if ismember(mode,[3,4]) 
               njump =1;
           elseif mode ==2;
               njump = n_randomly;
           end          
            id_greater = vsel > n_var;
            
            if ~isempty(vv(id_greater))
                   vsel(id_greater) = vsel(id_greater)-n_var;
            end
            
            if (mode ==4)
                    
                    irand = zeros(1,n_randomly-1);

                        for ii=1:n_randomly-1
                            isrepeat = false;
                            while ~isrepeat
                                irr =round(1+(n_var-1)*rand(1,1));  %generate N random number in the interval (a,b); where a=1, b=n_var(n_var=2) and N=1 with the formula r=a+(b-a).*rand(N,1)
                                if ~ismember(irr,[vsel(1) irand]);  
                                    isrepeat = true;
                                    irand(ii) = irr;
                                end
                            end
                        end

                    vsel(2:n_randomly) = irand;
            end %end if (mode==4)
            considered(vsel) = true;
            vsel = vsel+njump;
       elseif mode ==1
           vsel = randperm(n_var);
           vsel = vsel(1:n_randomly);
           considered(vsel) = true;
       end %end if ismember(mode,[2,3,4]) 
            
    end %end "VariableSelect1()"

    function x = h_function(x_bar,s1,s2,x_p)
        H=x_bar*(1.d0-exp(-x_p*s1))+(1.0d0-x_bar)*exp(-(1.d0-x_p)*s2);              
        H0=(1.d0-x_bar)*exp(-s2)  ;
        H1= x_bar*exp(-s1)  ;
        x=H+H1*x_p+H0*(x_p-1.d0);
    end
    function f_fit = static_penalty(fxin,gxin)
        ggxin = gxin;
        ggxin(gxin<0) = 0;
        f_fit = fxin+1e2*sum(ggxin.^2);
    end
   
end %function MVMO_old_nakawiro










