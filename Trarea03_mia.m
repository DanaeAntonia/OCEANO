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
%% calculamos tau 
rho_mar = 1025; %kg/m^3
rho_aire = 1.2; %kg/m^3
cd = 1.3e-3;
for i=1:3
tau_y(:,i)= rho_aire*cd*y_c(:,i).*mag_W(:,i);
end
%% Calculamos transporte
lat=[-37,-30,-21];
for i=1:3
M_x(:,i)= tau_y(:,i) ./ (rho_mar*(2*7.292e-5*sind(lat(i))));
end
%% Graficamos 
lat=[37,30,21];
fechas = datenum(fecha(:,1),fecha(:,2),fecha(:,3));
for i=1:3 
figure()
plot(fechas(1:365),M_x(1:365,i))
datetick('x','yyyy-mm')
ylabel('Transporte [m^2/s]')
xlabel('Timepo')
axis tight
grid on
title(['Serie de tiempo ', num2str(lat(i)),'°S']);
end
%% 2) Calcule el ciclo anual del transporte de Ekman promediando todos los valores de 
% enero, todos los de febreros y así hasta diciembre. Grafique el ciclo anual. Considere el 
% transporte de Ekman como un Índice de Surgencia (IS) y analice la figura de los ciclos 
% anuales en términos de la variabilidad anual de la surgencia y como estos varían en 
% las diferentes latitudes. 
%En la columna dos de las fechas tenemos los meses. 
% Extraer las filas para enero (mes 1)
for j=1:12
meses_37(:,j)= nanmean(M_x(fecha(:, 2) == j, 1));
meses_30(:,j)= nanmean(M_x(fecha(:, 2) == j, 2));
meses_21(:,j)= nanmean(M_x(fecha(:, 2) == j, 3));
end 
figure()
plot(meses_37,'g','LineWidth',2)
hold on 
plot(meses_30,'b','LineWidth',2)
plot(meses_21,'r','LineWidth',2)
title('Transporte de Ekman promediado por meses')
legend('37°S','30°S','21°S')
ylabel('Trasporte de Ekman [m^2/s]')
grid on 
axis tight
% Establecer los puntos para los meses en el eje x
xticks(1:length(meses_37))
% Agregar los nombres de los meses como etiquetas en el eje x
nombres_meses = {'Enero', 'Febrero', 'Marzo', 'Abril', 'Mayo', 'Junio', 'Julio', 'Agosto', 'Septiembre', 'Octubre', 'Noviembre', 'Diciembre'};
xticklabels(nombres_meses)
%% 3) Considere que el transporte hacia fuera de la costa tiene lugar en una capa superficial 
% de profundidad hE = 20 m y que este transporte es compensado por un transporte vertical 
% (hacia la capa superficial) en una distancia L desde la costa (ver esquema en figura 
% abajo y sus apuntes de clases). La distancia L puede ser relacionada con el radio de 
% deformación interno de Rossby dado por LR = ci / f , donde ci es la velocidad de fase de 
% una onda larga interna de gravedad. Use ci = 3 m s-1
% y f correspondiente a la latitud dada
% en cada caso.
lat=[-37,-30,-21];
for i=1:3 
LR=3/(2*7.292e-5*sind(lat(i)));
W(:,i)=(M_x(:,i)./(LR))*60*60*24;
clear LR
end 
% Definir los intervalos del histograma
intervalos = -5:1:10; % Rango de -5 m/día a 10 m/día con incremento de 1 m/día
% Obtener las fechas correspondientes a cada estación del año
verano_fechas =find( fecha(:,2) >= 1 & fecha(:,2) <= 3);
otono_fechas = find(fecha(:,2) >= 4 & fecha(:,2) <= 6);
invierno_fechas = find(fecha(:,2) >= 7 & fecha(:,2) <= 9);
primavera_fechas = find(fecha(:,2) >= 10 & fecha(:,2) <= 12);
% Calcular las velocidades verticales para cada estación del año y localidad
verano = W(verano_fechas, :);
otono = W(otono_fechas, :);
invierno = W(invierno_fechas, :);
primavera = W(primavera_fechas, :);
% Graficar histogramas para cada estación del año y localidad
lat=[37,30,21];
for i=1:3
figure()
subplot 221 
 histogram(verano(:, i), intervalos, 'Normalization', 'probability');
    title([num2str(lat(i)), '°S- Verano']);
    xlabel('Velocidad vertical (m/día)');
    ylabel('Porcentaje');
    ytix = get(gca, 'YTick');
    set(gca, 'YTick',ytix, 'YTickLabel',ytix*100);
subplot 222
 histogram(invierno(:, i), intervalos, 'Normalization', 'probability');
      title([num2str(lat(i)), '°S- Invierno']);
    xlabel('Velocidad vertical (m/día)');
    ylabel('Porcentaje');
    ytix = get(gca, 'YTick');
    set(gca, 'YTick',ytix, 'YTickLabel',ytix*100);
subplot 223
 histogram(otono(:, i), intervalos, 'Normalization', 'probability');
  title([num2str(lat(i)), '°S- Otoño']);
    xlabel('Velocidad vertical (m/día)');
    ylabel('Porcentaje');
    ytix = get(gca, 'YTick');
    set(gca, 'YTick',ytix, 'YTickLabel',ytix*100);
subplot 224
 histogram(primavera(:, i), intervalos, 'Normalization', 'probability');
   title([num2str(lat(i)), '°S- Primavera']);
    xlabel('Velocidad vertical (m/día)');
    ylabel('Porcentaje');
    ytix = get(gca, 'YTick');
    set(gca, 'YTick',ytix, 'YTickLabel',ytix*100);
    set(gcf,'color','w')
end
%
%% 4) Considere la serie de 37°S durante el verano del año 2000 (1 enero al 31 marzo)
dosmil=find(fecha(:,1)==2000 & fecha(:,2)>=1 & fecha(:,2)<=3 & fecha(:,3)>=1 &fecha(:,3)<=31);
%% a) Grafique (puede usar un gráfico de barras) la serie de tiempo del esfuerzo del viento
%paralelo a la costa (promedios diarios de Tau_y) durante el periodo.
tau=tau_y(dosmil,1);
% Crear gráfico de la serie de tiempo
dos= datenum(fecha(dosmil,1),fecha(dosmil,2),fecha(dosmil,3));
figure()
bar(dos, tau, 'FaceColor', 'b'); % Graficar la serie de tiempo como barras azules
xlabel('Fecha');
ylabel('Esfuerzo del Viento Paralelo a la Costa (Tau_y)');
title('Serie de Tiempo del Esfuerzo del Viento Paralelo a la Costa');
datetick('x','mm/dd')
axis tight
grid on
grid on; % Mostrar las líneas de la grilla
set(gca, 'FontName', 'Arial', 'FontSize', 12); % Personalizar la fuente y el tamaño del texto en los ejes
%% b) Calcule el impulso I del viento para los eventos del 7 al 11 de enero de 2000
fechas = datenum(fecha(:,1),fecha(:,2),fecha(:,3));
salo=find(fecha(:,1)==2000 & fecha(:,2)==1 & fecha(:,3)>=7 & fecha(:,3)<=11);
tau_I=tau_y(salo,1);
%los límites de tiempo inicial y final 'Di' y 'Df'
% Definir el intervalo de tiempo correspondiente
t_I = fechas(salo)*(3600*24); % Asegúrate de que el vector 't' contenga el tiempo correspondiente
% Calcular la integral definida utilizando el método del trapecio
impulso_E1 = (1/(50*1025))*trapz(t_I, tau_I);
disp(['El impulso del viento para los eventos del 7 al 11 de enero de 2000 es: ', num2str(impulso_E1)]);
%% Calcular el impulso I del viento para el evento del 13 al 27 de enero de 2000 
talv=find(fecha(:,1)==2000 & fecha(:,2)==1 & fecha(:,3)>=13 & fecha(:,3)<=27);
tau_I_E2=tau_y(talv,1);
%los límites de tiempo inicial y final 'Di' y 'Df'
% Definir el intervalo de tiempo correspondiente
t_I_E2 = fechas(talv)*3600*24; % Asegúrate de que el vector 't' contenga el tiempo correspondiente
% Calcular la integral definida utilizando el método del trapecio
impulso_E2 = (1/(50*1025))*trapz(t_I_E2, tau_I_E2);
disp(['El impulso del viento para los eventos del 13 al 27 de enero de 2000 es: ', num2str(impulso_E2)]);
%% Graficamos el perfil de velocidad
g=(2.5)^2 / 50;
A_E1=0.5476 * sqrt(50 / 0.1250);
LR=3/(2 * 7.292e-5 * sind(-37));
x=linspace(0, LR, 100);
R=(2.5)/(2 * 7.292e-5 * sind(-37));
for i = 1:length(x)
 v_E1(i)=A_E1*sqrt(0.1250/50)*exp(-x(i)/R);
end
d=3.0380/(2*7.292e-5*sind(-37))-R;
A_E2=50*exp(d/R);
for i = 1:length(x)
 v_E2(i)=A_E2*sqrt(0.1250 / 50)*exp(-x(i)/R);
end
figure()
plot(-1*flip(x), flip(v_E1),'b','LineWidth',3)
hold on
plot(-1*flip(x),flip(v_E2),'r','LineWidth',3)
grid on 
xlabel('Distancia desde la costa hacia afuera (metros)');
ylabel('Velocidad del chorro costero (m/s)');
title('Perfil del chorro costero');
legend('E1 7 al 11', 'E2 13 al 27')
axis tight
