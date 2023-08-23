%% A partir de los datos de velocidad del viento en el archivo VientosCosteros.mat
% calcule las series de tiempo del transporte de Ekman hacia fuera de la costa (dirección x)
clear all 
clc
load('VientosCosteros.mat'); 
alfa=zeros(3769,3);
for i=1:3
alfa(:,i)= atand(Ui(:,i)./Vi(:,i)); %Calculamos el angulo entre el vector de velocidad y la componente V de la velocidad del viento
aux=find(Ui(:,i)>0 & Vi(:,i)<0 );
alfa(aux)=alfa(aux)+180;
aux2=find(Ui(:,i)<=0 & Vi(:,i) < 0);
alfa(aux2,i)=alfa(aux2(i))-180;
clear aux aux2
mag_W(:,i) = ((Vi(:,i)).^2 + (Ui(:,i)).^2).^(1/2); %Magnitud del viento 
end 
%En alfa la column 1 es 37 la segunda 30 la tercera 21 
%Calculamos los nuevos marcos de referencia alineados con la costa 
y_c=zeros(3769,3);
x_c=zeros(3769,3);
y_c(:,1) = (mag_W(:,1)).*cosd(alfa(:,1)-30); %los 30 grados es el ángulo de la costa respecto al norte, esto se saca de Google Earth  
x_c(:,1) = mag_W(:,1).*sind(alfa(:,1)-30); 
y_c(:,2) = mag_W(:,2).*cosd(alfa(:,2)-5);%los 5 sacado de Google Earth
x_c(:,2) = mag_W(:,2).*sind(alfa(:,2)-5);
y_c(:,3) = mag_W(:,3).*cosd(alfa(:,3)-3);%los 3 los sacas de Google Earth 
x_c(:,3) = mag_W(:,3).*sind(alfa(:,3)-3);