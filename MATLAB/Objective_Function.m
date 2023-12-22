function [f,g,xn_out] = Objective_Function(xx_yy1,iss)

global parameter
global p_des p_act P_base eff_ramp option raw_ACE t const1

% ----------------  Denormalize to the original range ---------------------
 Scaling = (parameter.x_max-parameter.x_min);    
 xx_yy = parameter.x_min+Scaling.* xx_yy1;
% -------------------------------------------------------------------------

%1: tracking logic
%2: rate limiting logic
%3: ACE logic

switch iss
    case 1 %tracking logic
        DELT=4;
        acum_ERR=0;
        acum_SUM2=0;
        size=length(t);
        SUM1=zeros(1,size);
        SUM2=zeros(1,size);
        ERR=zeros(1,size);
        for i=1:size
            if i==1
                SUM1(i)=p_des(i);
                SUM2(i)=0;
                ERR(i)=p_act(i)-SUM1(i);
                acum_ERR=abs(abs(ERR(i))+acum_ERR);
                acum_SUM2=abs(abs(SUM2(i))+acum_SUM2);
            else
                SUM1(i)=SUM1(i-1)+(DELT/xx_yy(1))*(p_des(i)-SUM1(i-1));
                ERR(i)=p_act(i)-SUM1(i);
                acum_ERR=abs(abs(ERR(i))+acum_ERR);
                SUM2(i)=SUM2(i-1)+(DELT/xx_yy(2))*(ERR(i)-SUM2(i-1));
                acum_SUM2=abs(abs(SUM2(i))+acum_SUM2);
            end
        end
       
       %Objective function: Minimize accumulated ERR
        f=acum_ERR;
        
        %Constraints
        switch option
            case 1 %UP
                r1=0;
                for i=1:size
                    if round(SUM1(i))>=P_base
                        break
                    else
                        if SUM1(i)>=p_act(i)
                            r1=r1+0;
                        else
                            r1=r1+1;
                        end
                    end
                end
                
                if r1>0
                    g=1;
                else
                    g=0;
                end
                
            case 2 %DOWN
                r1=0;
                for i=1:size
                    if round(SUM1(i))<=P_base
                        break
                    else
                        if SUM1(i)<=p_act(i)
                            r1=r1+0;
                        else
                            r1=r1+1;
                        end
                    end
                end
                
                if r1>0
                    g=1;
                else
                    g=0;
                end    
        end
        
    case 2 %rate limiting logic
        DELT=4;
        size=length(t);
        SUM3=zeros(1,size);
        SUMP=zeros(1,size);
        SUM4=zeros(1,size); 
        for i=1:size
            if i==1
                SUM3(i)=p_act(i);
                SUMP(i)=0;
                SUM4(i)=SUMP(i);             
            else
                SUM3(i)=SUM3(i-1)+(DELT/xx_yy(1))*(p_act(i)-SUM3(i-1));
                SUMP(i)=(p_act(i)-SUM3(i-1))*60/xx_yy(1);
                SUM4(i)=SUM4(i-1)+(DELT/xx_yy(2))*(SUMP(i)-SUM4(i-1));
            end
        end
        
        %Objective function
        ERR=abs(eff_ramp-max(abs(SUMP)));
        f=ERR;
        %Constraint
        g=0;
            
    case 3 %ACE logic
        %Moving average calculation
        win_size=xx_yy(1); %window size, number of samples 
        p=length(raw_ACE);
        
        for l=1:p
            if l<win_size
                n=l;
            else
                n=win_size;
            end
            for i=1:n
                aux(i)=raw_ACE(l-i+1);
            end
            ma(l)=(1/n)*sum(aux);
        end
        
        %Filtered ACE calculation
        % exponential ACE
        size=length(raw_ACE);
        ACE_exp=zeros(1,size);
        DELT=4;
        XWMATC=xx_yy(2); 
        for i=1:length(raw_ACE)
            if i==1
                ACE_exp(i)=raw_ACE(i);
            else
                ACE_exp(i)=ACE_exp(i-1)+(DELT/XWMATC)*(raw_ACE(i)-ACE_exp(i-1));
            end
        end
        %ACE moving average
        ACE_mov_avg=ma;
        
        %Filtered ACE
        G_avg=0.5; 
        G_exp=0.5; 
        FACE=G_avg*ACE_mov_avg+G_exp*ACE_exp;
        
        %Accuracy
        for i=1:p
            aux2(i)=abs(raw_ACE(i)-FACE(i));
        end
        acc=(1/p)*sum(aux2);
        
        %Smoothness
        for i=3:p
            aux3(i)=abs(FACE(i)-2*FACE(i-1)+FACE(i-2));
        end
        smo=(1/p)*sum(aux3);
        
        %Objective function
        ERR=acc;
        f=ERR;
           
        %Constraints
        r1=abs(smo-const1);
        if r1<0.01
            g=0;
        else
            g=1;
        end       
 end
% ------------------------  Normalize again -------------------------------
xn_out = (xx_yy-parameter.x_min)./Scaling;
% -------------------------------------------------------------------------

end



