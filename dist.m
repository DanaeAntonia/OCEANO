function d= dist(lat1,lat2,lon1,lon2) 
%defino radio de la tierra 6371*1e3 m 
Rt=6371*1e3; 
%transformo grados a radiane
lat1_rad = deg2rad(lat1); 
lat2_rad = deg2rad(lat2);
lon1_rad = deg2rad(lon1);
lon2_rad = deg2rad(lon2);
% Calcular la distancia utilizando la fórmula del coseno del haversine
d= 2 * Rt * asin(sqrt(sin((lat2_rad - lat1_rad) / 2).^2 + cos(lat1_rad) * cos(lat2_rad) * sin((lon2_rad - lon1_rad) / 2).^2));
end 


