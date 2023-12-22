close all
clear ERR_mod
global p_des p_act P_base option t

DELT=4;

%Ingreso de parámetros
T1=9;  %Valor óptimo
T2=60; %Valor óptimo entre: 40-80 [s]

acum_ERR=0;
acum_SUM2=0;
size=length(t);
SUM1=zeros(1,size);
SUM2=zeros(1,size);
ERR=zeros(1,size);
for i=1:size
    if i==1
        SUM1(i)=p_des(i);
        ERR(i)=0;
        SUM2(i)=0;
    else
        SUM1(i)=SUM1(i-1)+(DELT/T1)*(p_des(i)-SUM1(i-1));
        ERR(i)=p_act(i)-SUM1(i-1);
        SUM2(i)=SUM2(i-1)+(DELT/T2)*(ERR(i)-SUM2(i-1));
    end
end

%Determinación K6 (Umbral de no seguimiento)
Percent=50;  %Porcentaje por encima del pico de SUM2 %Mínimo 10%
K6_aux=round(max(abs(SUM2))*(1+Percent/100));

switch option
    case 1 %UP
        K6=-K6_aux;
    case 2 %DOWN
        K6=K6_aux;
end

K6_plot=0*t+K6;

%Determinación KD (Banda muerta para ignorar el error de seguimiento)
switch option
    case 1 %UP
        for i=1:size
            if round(p_act(i))>=P_base
                break
            else
                ERR_mod(i)=ERR(i);
            end
        end
    case 2 %DOWN
        for i=1:size
            if round(p_act(i))<=P_base
                break
            else
                ERR_mod(i)=ERR(i);
            end
        end
end
Percent=5; %Porcentaje por debajo de la diferencia más pequeña entre PWGENX y SUM1
KD=floor(min(abs(ERR_mod))*(1-Percent/100));

%Print and Plots
fprintf('\nK6 = %.4f\n',K6)
fprintf('\nKD = %.4f\n',KD)
figure
plot (t,p_act,t,p_des,t,SUM1,t,SUM2,t,K6_plot,'LineWidth',2)
legend('Pot. actual','Pot.deseada','SUM1','SUM2','K6')
figure
plot (t,SUM2,t,K6_plot,'LineWidth',2)
legend('SUM2','K6')