% Valores
f= abs(2*7.292e-5*sind(-37));
Av = 1e-2; %[m^2/s]
r0 = 1025;
d= sqrt(2*Av/f);
z= -(0:0.5:130)';
T= 0.1; %[Pa]
% Velocidades de Ekman
u = T/(r0*sqrt(f*Av))*exp(z/d).*cos(-z/d +pi/4);
v = T/(r0*sqrt(f*Av))*exp(z/d).*sin(-z/d +pi/4);

%% defino vectores desde el centro de la espiral
x = randn(size(u)); %defino vectores
y = randn(size(v));
center_x = zeros(size(x)); %defino desde el centro
center_y = zeros(size(y));
end_x = center_x + u;
end_y = center_y + v;
idx = 1:5:261; % definir solo algunos vectores
%% GRAFICO
figure()
subplot(121)
quiver(center_x(idx),center_y(idx),u(idx),v(idx),'k');
xlabel('u [m/s]')
ylabel('v [m/s]')
title('Espiral de Ekman')
grid on
set(gcf,'color','w')
axis([-0.08 0.08 -0.08 0.08])
subplot(122)
plot(u,z,'b','LineWidth',2)
hold on
plot(v,z,'r','LineWidth',2)
legend('u','v')
xlabel('velocidades u,v [m/s]')
ylabel('Profundidad [m]')
title('Perfil de velocidades')
grid on
set(gcf,'color','w')

%% 
figure(2) %grafico
quiver(center_x(idx),center_y(idx),u(idx),v(idx),'k');
xlabel('u [m/s]')
ylabel('v [m/s]')
title('Espiral de Ekman')
grid on
set(gcf,'color','w')
%% realizamos perfil de velocidad 
figure(3)
plot(u,z,'b','LineWidth',2)
hold on
plot(v,z,'r','LineWidth',2)
legend('u','v')
xlabel('velocidades u,v [m/s]')
ylabel('Profundidad [m]')
title('Perfil de velocidades')
grid on
set(gcf,'color','w')
%% meto todo en un subplot

%% cambio los parametros para ver que apsa
Avv=[1e-2 ; 2e-1 ; 1 ; 0 ; 10]
r00=[1025 ; 500; 250; 1500; 50]
dd= sqrt(2*Avv/f);
TT=[0.1 ; 0.5; 1; 10; 50]
% Velocidades de Ekman
for i=1:5
   uu = T/(r0*sqrt(f*Avv(i)))*exp(z/d).*cos(-z/d +pi/4);
   vv = T/(r0*sqrt(f*Avv(i)))*exp(z/d).*sin(-z/d +pi/4);
   x = randn(size(uu)); %defino vectores
   y = randn(size(vv));
   center_x = zeros(size(x)); %defino desde el centro
   center_y = zeros(size(y));
   idx = 1:5:261;
   figure(i)
   subplot(121)
    quiver(center_x(idx),center_y(idx),uu(idx),vv(idx),'k');
   xlabel('u [m/s]')
   ylabel('v [m/s]')
    title('Espiral de Ekman')
    grid on
    set(gcf,'color','w')
    axis equal
    subplot(122)
    plot(uu,z,'b','LineWidth',2)
    hold on
    plot(vv,z,'r','LineWidth',2)
    legend('u','v')
    xlabel('velocidades u,v [m/s]')
    ylabel('Profundidad [m]')
    title('Perfil de velocidades')
    grid on
    set(gcf,'color','w')
end

