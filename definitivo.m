clear all clc
%% Grafique los campos originales del rotor para las distintas estaciones 
%del año.
load('WSCurlPSO_12.mat'); %cargo el rotor del viento
%separa por estaciones
verano=WSCurl(:,:,1:3); 
otono=WSCurl(:,:,4:6); 
invierno=WSCurl(:,:,7:9); 
primavera=WSCurl(:,:,10:12);
%calculo promedio para obtener el rotor para cada estación 
%multiplico por 10^-7 ya que el rotor esta con ese orden.
ver=mean(verano,3); 
oto=mean(otono,3);
inv=mean(invierno,3);
pri=mean(primavera,3);
% Definir los colores base
negro = [0, 0, 0]; % Negro
morado_oscuro = [0.5, 0, 0.5]; % Morado oscuro
morado_claro = [0.8, 0.4, 0.8]; % Morado claro
blanco = [1, 1, 1]; % Blanco
turquesa_claro = [0.4, 0.8, 0.8]; % Turquesa claro
azul_marino = [0, 0, 0.5]; % Azul marino

% Definir los puntos de referencia y los colores correspondientes
x = [0, 0.2, 0.4, 0.6, 0.8, 1]; % Puntos de referencia
colores = [azul_marino; turquesa_claro; blanco; morado_claro; morado_oscuro; negro]; % Colores correspondientes

% Crear una paleta de colores interpolada
num_colores = 8; % Número de colores en la paleta
mapa_colores = interp1(x, colores, linspace(0, 1, num_colores));


[lati,long] = meshgrid(double(lat1),double(lon1)); %grilla para graficar
figure()
subplot 221
m_proj('Robinson','lon',[min(long(:)) max(long(:))],'lat',[min(lati(:)) max(lati(:))])
m_pcolor(lon1, lat1, ver')%, 'LineWidth', 0.1);
m_coast('line', 'color', 'k', 'linewidth', 2);
m_grid('box', 'fancy', 'tickdir', 'in');
colormap(mapa_colores);
colorbar;
title('Verano');
set(gca, 'FontSize', 12);
caxis([-4 4])
hold on
m_contour(lon1,lat1,ver',[1 1],'LineColor', 'k','linewidth',0.5)
m_contour(lon1,lat1,ver',[0 0],'LineColor', 'k','linewidth',0.5)
m_contour(lon1,lat1,ver',[2 2],'LineColor', 'k','linewidth',0.5)
m_contour(lon1,lat1,inv',[-1 -1],'LineColor', 'k','linewidth',0.5)
subplot 222
m_proj('Robinson','lon',[min(long(:)) max(long(:))],'lat',[min(lati(:)) max(lati(:))])
m_pcolor(lon1, lat1, oto')%, 'LineWidth', 0.1);
m_coast('line', 'color', 'k', 'linewidth', 2);
m_grid('box', 'fancy', 'tickdir', 'in');
colormap(mapa_colores);
colorbar;
title('Otoño');
caxis([-4 4])
set(gca, 'FontSize', 12);
hold on
m_contour(lon1,lat1,oto',[1 1],'LineColor', 'k','linewidth',0.5)
m_contour(lon1,lat1,oto',[0 0],'LineColor', 'k','linewidth',0.5)
m_contour(lon1,lat1,oto',[2 2],'LineColor', 'k','linewidth',0.5)
m_contour(lon1,lat1,inv',[-1 -1],'LineColor', 'k','linewidth',0.5)
subplot 223
m_proj('Robinson','lon',[min(long(:)) max(long(:))],'lat',[min(lati(:)) max(lati(:))])
m_pcolor(lon1, lat1, pri')
m_coast('line', 'color', 'k', 'linewidth', 2);
m_grid('box', 'fancy', 'tickdir', 'in');
colormap(mapa_colores);
colorbar;
title('Primavera');
set(gca, 'FontSize', 12);
caxis([-4 4])
hold on
m_contour(lon1,lat1,pri',[1 1],'LineColor', 'k','linewidth',0.5)
m_contour(lon1,lat1,pri',[0 0],'LineColor', 'k','linewidth',0.5)
m_contour(lon1,lat1,pri',[2 2],'LineColor', 'k','linewidth',0.5)
m_contour(lon1,lat1,inv',[-1 -1],'LineColor', 'k','linewidth',0.5)
subplot 224
set(gca, 'FontSize', 12);
m_proj('Robinson','lon',[min(long(:)) max(long(:))],'lat',[min(lati(:)) max(lati(:))])
m_pcolor(lon1, lat1, inv')
m_coast('line', 'color', 'k', 'linewidth', 2);
m_grid('box', 'fancy', 'tickdir', 'in');
colormap(mapa_colores);
colorbar;
title('Invierno');
set(gca, 'FontSize', 12);
set(gcf,'color','w')
caxis([-4 4])
hold on
m_contour(lon1,lat1,inv',[1 1],'LineColor', 'k','linewidth',0.5)
m_contour(lon1,lat1,inv',[0 0],'LineColor', 'k','linewidth',0.5)
m_contour(lon1,lat1,inv',[2 2],'LineColor', 'k','linewidth',0.5)
m_contour(lon1,lat1,inv',[-1 -1],'LineColor', 'k','linewidth',0.5)
sgtitle('Rotor del esfuerzo del viento [Pa/m por 10^{-7}] para cada estación del año','FontSize',16)
%% Grafique el campo de wE en todo el dominio (excluya de sus cálculos de w la banda
%ecuatorial entre -3°S 3°N). Destaque en la figura el contorno wE = 0
rho=1025; %kg/m^3
a1=lat1<3 & lat1>-3;
WSCurl(:,a1,:)=NaN; %le saco la banda de -3 a 3
rotor=mean(WSCurl,3); %Calculo promedio anual del rotor
f=2*7.29*10^(-5).*sind(lat1);% calculo el parametro de coriolis para las diferentes latitudes
w=rotor'./(rho.*f)*10^(-7); % calculo las velocidades verticales %m/s
%% Grafico
figure()
m_proj('Robinson','lon',[min(long(:)) max(long(:))],'lat',[min(lati(:)) max(lati(:))])
m_pcolor(lon1, lat1, w); shading flat;
colorbar
colormap('jet')
caxis([-0.3e-5 0.3e-5])
hold on 
m_contour(lon1,lat1,w,[0 0],'LineColor', 'k','linewidth',2)
m_coast('line', 'color', 'k', 'linewidth', 1);
m_grid('box', 'fancy', 'tickdir', 'in');
title('Velocidades verticales asociadas al rotor del viento [m/s]');
set(gca, 'FontSize', 12);
set(gcf,'color','w')
%% Calcule la componente meridional del transporte de Sverdrup (V ) en todo el dominio
%usando los datos promedios de todo el periodo y la relación vista 
clear all clc
load('WSCurlPSO_12.mat')
Curl=mean(WSCurl,3)*1e-7; % promedio de rotor del viento
B= 2*1e-11;
V= (Curl' ./ (1025 .* B));
%% Graficar
[lati,long] = meshgrid(double(lat1),double(lon1)); %grilla para graficar
figure()
m_proj('Robinson','lon',[min(long(:)) max(long(:))],'lat',[min(lati(:)) max(lati(:))])
m_pcolor(lon1, lat1, V); shading flat;
colorbar
colormap('jet')
caxis([-10 10]) 
hold on 
m_contour(lon1,lat1,V,[0 0],'LineColor', 'k','linewidth',2)
m_coast('line', 'color', 'k', 'linewidth', 1);
m_grid('box', 'fancy', 'tickdir', 'in');
title('Transporte meridional [m^2/s]');
set(gca, 'FontSize', 12);
set(gcf,'color','w')


%% integramos el transporte de svedrup a lo largo de longitudes 
%para una latitud fija
clear all clc
load('WSCurlPSO_12.mat');
Curl=mean(WSCurl,3)*1e-7; % promedio de rotor del viento
y0 = [-20.1250000 -30.1250000 -40.1250000]; % Latitudes fijas
%deberia dejarlo constante 
B=2*1e-11; 
%B = 2 * 7.29 * 10^(-5) * cosd(y0) / (6*10^(6)); % Cálculo de B en función de y0
% Definir el rango de longitud en el que deseas integrar
lon22= flip(lon1); % valores de lon1 desde costa hacia afuera
%debo converitir las longitudes en distancias, para eso utilizaremos el
%metodo del coseno del haversine, con la funcion dist, creada por mi 
for i=1:length(lat1)
    c=diff(lat1);
  for   j=1:length(lon1)
      d=diff(lon1);
      lat2=lat1(i)+c(1);
      lon2=lon1(j)+d(1);
      lat=lat1(i);
      lon=lon1(j);
      distancia(j,i)=dist(lat,lat2, lon,lon2);
  end 
end 
% Inicializar el vector para almacenar los resultados de Vt
Vt = zeros(size(lon1));
% Calcular Vt para cada punto (x(i), y0)
for j=1:3
for i = 1:length(lon1)
    Vt(i,j) = Curl(i, find(lat1 == y0(j), 1)) / (1025 * B);
end
end
v2=flipud(Vt); % cambio las filas, desde la costa hacia afuera 
indices_nan = isnan(v2); %encuentro posiciones NaN 
v2(indices_nan) = 0;% Reemplazar los elementos NaN por 0
%% para 20°S 
lon_20=distancia(1,find(lat1==y0(1)));
integral_20=cumtrapz(lon22, v2(:,1))*lon_20(1);
%% para 30°S
lon_30=distancia(1,find(lat1==y0(2)));
integral_30=cumtrapz(lon22, v2(:,2))*lon_30(1);
%% para 40°S
lon_40=distancia(1,find(lat1==y0(3)));
integral_40=cumtrapz(lon22, v2(:,3))*lon_40(1);
%% generalizo todo a una matriz 
integral=[integral_20 , integral_30,  integral_40];
%% grafico
figure()
plot(lon1,-1*flipud(integral(:,1))*1e-6,'b','LineWidth',3)
hold on
plot(lon1,-1*flipud(integral(:,2))*1e-6,'r','LineWidth',3)
plot(lon1,-1*flipud(integral(:,3))*1e-6,'k','Linewidth',3)
line([min(lon1), max(lon1)], [0, 0], 'LineStyle', '--', 'Color', 'g')
grid on 
xlabel('Longitud')
ylabel('Transporte [Sv]')
title('Transporte de svedrup integrado en longitud')
legend('lat 20°S' , 'lat 30°S', 'lat 40°S')
set(gcf,'color','w')
%% calcular el transporte zonal U(x, y) 
clear all clc 
load('WSCurlPSO_12.mat');
Curl=mean(WSCurl ,3); 
B= 2*1e-11;
V= Curl' ./ (1025 .* B) .*10^(-7);
% Calcular dv/dy iterativamente
dv_dy = zeros(size(V));
VV=flipud(V);
lat2=flip(lat1);
%voy a pasar a distancias mis latitudes para tener consistencia 
% dis_lat=diff(lat1)*111*1e3; %asumimos que un grado corresponde a 111 km aprox
% dis_lat(240)=dis_lat(1);%solo para tener las mismas dimensiones
%cambio NaN por ceros.
% indices_nan = isnan(VV);
% Reemplazar los elementos NaN por 0
% VV(indices_nan) = 0;
for i = 1:size(V, 2)
    dv_dy(:, i) = gradient(VV(:, i), lat2);%*dis_lat(1);
end
%% ahora integro 
lon2=flip(lon1); %doy vuelta las longitudes 
U= cumtrapz(lon2, dv_dy')';
%debo multiplicar por la longitud en distancia para cada latitud. 
for i=1:length(lat1)
    c=diff(lat1);
  for   j=1:length(lon1)
      d=diff(lon1);
      lat2=lat1(i)+c(1);
      lon2=lon1(j)+d(1);
      lat=lat1(i);
      lon=lon1(j);
      distancia(j,i)=dist(lat,lat2, lon,lon2);
  end 
end 
distancia=distancia';
for i=1:200
    UU(:,i)=U(:,i)*distancia(1,i);
end
%% Grafico
% Definir los colores base
negro = [0, 0, 0]; % Negro
rosado_intenso = [1, 0, 0.5]; % Rosado intenso
rosado_clarito = [1, 0.6, 0.8]; % Rosado clarito
blanco = [1, 1, 1]; % Blanco
celeste = [0.6, 0.8, 1]; % Celeste
azul_intenso = [0, 0, 0.5]; % Azul intenso

% Definir los puntos de referencia y los colores correspondientes
x = [0, 0.2, 0.4, 0.6, 0.8, 1]; % Puntos de referencia
colores = [azul_intenso; celeste; blanco; rosado_clarito; rosado_intenso; negro]; % Colores correspondientes

% Crear una paleta de colores interpolada
num_colores = 19; % Número de colores en la paleta
mapa_colores = interp1(x, colores, linspace(0, 1, num_colores));

[lati,long] = meshgrid(double(lat1),double(lon1)); %grilla para graficar
figure()
m_proj('Robinson','lon',[min(long(:)) max(long(:))],'lat',[min(lati(:)) max(lati(:))])
m_pcolor(lon1, lat1, flipud(UU)*1e-5); shading flat;
colorbar
colormap(mapa_colores)
caxis([-15 15]) 
hold on 
m_contour(lon1,lat1,flipud(UU),[0 0],'LineColor', 'k','linewidth',1)
m_coast('line', 'color', 'k', 'linewidth', 1);
m_grid('box', 'fancy', 'tickdir', 'in');
title('Transporte zonal [Sv]');
set(gca, 'FontSize', 12);
set(gcf,'color','w')

%% Calcular la función Corriente
clear all clc 
load('WSCurlPSO_12.mat');
Curl=mean(WSCurl ,3); 
B= 2*1e-11;
V= Curl' ./ (1025 .* B) .*10^(-7); %calculo transporte meridional
VV=flipud(V');
lat11=lat1;%invierto para tener desde la costa hacia afuera.
%debo calcular la funcion para una latitud fija. 
lon22=flip(lon1);%invierto para tener desde la costa hacia afuera
%transformo lon en distancias 
for i=1:length(lat1)
    c=diff(lat1);
  for   j=1:length(lon1)
      d=diff(lon1);
      lat2=lat1(i)+c(1);
      lon2=lon1(j)+d(1);
      lat=lat1(i);
      lon=lon1(j);
      distancia(j,i)=dist(lat,lat2, lon,lon2);
  end 
end 
indices_nan = isnan(VV);
% Reemplazar los elementos NaN por 0
VV(indices_nan) = 0;
for i=1:length(lat11)
    phi(:,i)=flipud((cumtrapz(distancia(1,i),VV(:,i)')));
end 
phi(indices_nan)=NaN;


%% Grafico phi la función corriente 
[lati,long] = meshgrid(double(lat1),double(lon1)); %grilla para graficar
figure()
m_proj('Robinson','lon',[min(long(:)) max(long(:))],'lat',[min(lati(:)) max(lati(:))])
m_pcolor(flipud(long),lati,phi*1e-6)
colormap('jet')
colorbar
hold on
caxis([-10 15])
for i=-10:5:15
m_contour(flip(lon1),lat1,(phi*1e-6)',[i i],'LineColor', 'k','linewidth',0.5);
h=m_contour (flip(lon1),lat1,(phi*1e-6)',[i i],'LineColor', 'k','linewidth',0.5);
clabel(h,'FontSize',10,'LabelSpacing',300);
end
%m_contour(flip(lon1),lat1,(phi*1e-6)',[0 0],'LineColor', 'k','linewidth',0.5)
m_coast('line', 'color', 'k', 'linewidth', 1);
m_grid('box', 'fancy', 'tickdir', 'in');
title('Función corriente transporte de Svedrup [Sv]');
set(gca, 'FontSize', 12);
set(gcf,'color','w')
