close all
global p_act eff_ramp t

DELT=4;

%Ingreso de parámetros
T3=33;  %Valor óptimo
T4=20;  %Valor óptimo entre 10-40 [s] 

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
        SUM3(i)=SUM3(i-1)+(DELT/T3)*(p_act(i)-SUM3(i-1));
        SUMP(i)=(p_act(i)-SUM3(i-1))*(60/T3);
        SUM4(i)=SUM4(i-1)+(DELT/T4)*(SUMP(i)-SUM4(i-1));
    end
end

%Determinación de la rampa
switch option
    case 1
        syms ramp(t);
        eff_ramp_s=eff_ramp/60; %[Mw/seg]
        ramp(t)=eff_ramp_s*t;
        t=zeros(1,size);
        for i=2:size
            t(i)=t(i-1)+4;
        end
        ramp=double(ramp(t));
        ramp_mod=0;
        for i=1:size
            if ramp(i)>=eff_ramp
                ramp_mod(i)=eff_ramp;
            else
                ramp_mod(i)=ramp(i);
            end
        end 
    case 2
        syms ramp(t);
        eff_ramp_s=eff_ramp/60; %[Mw/seg]
        ramp(t)=-eff_ramp_s*t;
        t=zeros(1,size);
        for i=2:size
            t(i)=t(i-1)+4;
        end
        ramp=double(ramp(t));
        ramp_mod=0;
        for i=1:size
            if ramp(i)<=-eff_ramp
                ramp_mod(i)=-eff_ramp;
            else
                ramp_mod(i)=ramp(i);
            end
        end 
end

%Determinación tiempo en el cual SUM4 debe alcanzar pico
aux2=max(abs(SUMP));
index=find(abs(SUMP)==aux2);
t_to_T4=t(index);

%Determinación KRATLM (Límite de velocidad a largo plazo o limitación de rampa)
Percent=10;  %Porcentaje por encima del pico de SUM4 %mínimo 5%
aux_KRATLM=max(abs(SUM4))*(1+Percent/100);
KRATLM=round(max(abs(SUM4))*(1+Percent/100));
switch option
    case 1
        KRATLM_plot=0*t+KRATLM;
    case 2
        KRATLM_plot=0*t-KRATLM;
end

%Print and Plots
fprintf('\nKRATLM = %.4f\n',KRATLM)
plot (t,p_act,t,SUM3,t,SUMP,t,SUM4,t,ramp_mod,'--',t,KRATLM_plot,'LineWidth',2)
legend('Pot. actual','SUM3','SUMP','SUM4','RAMPA','KRATLM')