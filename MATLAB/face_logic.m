clear;
global t raw_ACE

%Ingreso de parámetros
win_size=36; %tamaño de la ventana, número de muestras 
XWMATC=70; 

%Cálculo moving average
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

%Cálculo ACE filtrado
% ACE exponencial
size=length(raw_ACE);
ACE_exp=zeros(1,size);
DELT=4;
 for i=1:length(raw_ACE)
    if i==1
        ACE_exp(i)=raw_ACE(i);
    else
        ACE_exp(i)=ACE_exp(i-1)+(DELT/XWMATC)*(raw_ACE(i)-ACE_exp(i-1));
    end
 end 
%Promedio móvil ACE
ACE_mov_avg=ma;

%ACE filtrado
G_avg=0.5; 
G_exp=0.5; 
FACE=G_avg*ACE_mov_avg+G_exp*ACE_exp;

%Accuracy
for i=1:p
    aux2(i)=abs(raw_ACE(i)-FACE(i));
end
acc=(1/p)*sum(aux2);
sumaux2=sum(aux2);
ERR=aux2-4;
sumERR=sum(abs(ERR));
Y=(ERR).^2;
f=trapz(t,Y);

%Smoothness
for i=3:p
    aux3(i)=abs(FACE(i)-2*FACE(i-1)+FACE(i-2));
end
smo=(1/p)*sum(aux3);

%MIN ACE
MIN_ACE_p=6+0*t;
MIN_ACE_n=-6+0*t;

%Print
fprintf('\nAcc = %.4f\n',acc)
fprintf('\nSmo = %.4f\n',smo)

%Figuras
figure
p=plot (t,raw_ACE,t,FACE,t,MIN_ACE_p,t,MIN_ACE_n,'LineWidth',2);
legend('raw ACE','FACE','MIN ACE +','MIN ACE -')
p(1).LineWidth = 1;
p(2).LineWidth = 2;
p(3).LineWidth = 0.5;
p(4).LineWidth = 0.5;