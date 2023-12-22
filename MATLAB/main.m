close all
clear xlRange xlRange2 xlRange3
global p_des p_act t P_base eff_ramp option raw_ACE const1

%1: tracking logic and rate limiting logic / parámetros originales rampa subida - Molino  
%2: tracking logic and rate limiting logic / parámetros originales rampa bajada - Molino  
%3: ACE filtering logic / parámetros originales FACE  
 
 option=3;
 switch option
     case 1 
         %Ingreso de parámetros
         eff_ramp=50; %[MW/min]
         
         %lectura datos Excel
         filename = 'parametros_originales_Molino.xlsx';
         sheet = 'Hoja1';
         xlRange = 'K76:K243';
         xlRange2 = 'C76:C243';
         xlRange3 = 'G76';
         p_des = xlsread(filename,sheet,xlRange);
         p_act = xlsread(filename,sheet,xlRange2);
         P_base= xlsread(filename,sheet,xlRange3);
         size=length(p_act);
         t=zeros(1,size);
         for i=2:size
             t(i)=t(i-1)+4;
         end
         plot (t,p_act,t,p_des,'LineWidth',2)
         legend('Pot. real','Pot. deseada')
     case 2
         %Ingreso de parámetros
         eff_ramp=50; %[MW/min]
         
         %lectura datos Excel
         filename = 'parametros_originales_Molino.xlsx';
         sheet = 'Hoja1';
         xlRange = 'K244:K454';
         xlRange2 = 'C244:C454';
         xlRange3 = 'G244';
         p_des = xlsread(filename,sheet,xlRange);
         p_act = xlsread(filename,sheet,xlRange2);
         P_base= xlsread(filename,sheet,xlRange3);
         size=length(p_act);
         t=zeros(1,size);
         for i=2:size
             t(i)=t(i-1)+4;
         end
         plot (t,p_act,t,p_des,'LineWidth',2)
         legend('Pot. real','Pot. deseada')
     case 3
         %Ingreso de parámetros
         const1=0.1; %Valor de suavizado constante para la restricción
         
         %Lectura datos Excel
         filename = 'parametros_originales_FACE.xlsx';
         sheet = 'Hoja1';
         xlRange = 'D2:D902';  %1 hora
         raw_ACE=xlsread(filename,sheet,xlRange);
         
         %Cálculo vector tiempo
         size=length(raw_ACE);
         t=zeros(1,size);
         for i=1:length(raw_ACE)
             if i==1
                 t(i)=0;
             else
                 t(i)=t(i-1)+4;
             end
         end
 end