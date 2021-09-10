clear all
close all

%% cargar info de la falla
% importante las coordenadas de la falla y observaciones en UTM
% 

% !cat  ../malla_slab/slab_geo.inp | awk -F ',' '{if($3>10000 && $3 ~ /[4-9]/) print $0}'| sed 's/,/ /g' > nodos_new.txt
% !cat  ../malla_slab/slab_geo.inp | awk -F',' '{if(NF>3 && $3<10000  && $2 ~ /[1-9]/ && $3>0.1) print $0}'| sed 's/,/ /g'> tri_new.txt


% tri=load('tri_new.txt');
% coord=load('nodos_new.txt');

tri=load('tri.txt');
coord=load('nodos.txt');
triangulos=tri;
nodos_triangulos =coord;

tri(:,1)=[];
xyz=coord(:,2:end);
n=size(tri,1); %numero de triangulos

% figura %triangulos falla 
figure (1)
trimesh(tri, xyz(:,1), xyz(:,2), xyz(:,3));
 
% coordenadas centros triangunlos % necesita funcion center.m
for  i=1:n
v1=[coord(tri(i,1),2), coord(tri(i,1),3), coord(tri(i,1),4) ];
v2=[coord(tri(i,2),2), coord(tri(i,2),3), coord(tri(i,2),4) ];
v3=[coord(tri(i,3),2), coord(tri(i,3),3), coord(tri(i,3),4) ];

centro(i,:)=center(v1, v2, v3, 'barycenter');
area(i,1)=triangleArea3d(v1, v2, v3);
end

%% estima distribucion de slip desde mapa o desde un archivo
slip_conocido=1;  %para seleccionar area de slip desde mapa
if slip_conocido==1

%% localizar donde ocurre un deslizamiento % necesita funcion selectdata.m
figure(1) 
plot(centro(:,1),centro(:,2),'x');
[pointslist,xselect,yselect] = selectdata;
close (figure(1))

% definir slip
%maximo en el centro disminuye con la distancia al centro

puntos=[centro(pointslist,:)];
pto_c=[mean(puntos(:,1)), mean(puntos(:,2)), mean(puntos(:,3))];

slip_max=1; %en metros 

%ponderar slip por distancias
for i=1:length(puntos)
    X=[ puntos(i,1), puntos(i,2) , puntos(i,3) ; pto_c  ];
    d(i)=pdist(X,'euclidean');
end

coefficients = polyfit([0, max(d)], [slip_max, 0], 1);
a = coefficients (1);
b = coefficients (2);

slip=a*d+b;

slipf=zeros(size(centro,1),1);
slipf(pointslist)=slip;   %slipf mismo size q numero de triangulos

save slipf.txt slipf -ASCII

else 
     slipf=load('slipf.txt');     
end


%% selecionar slip por profundidad
%  cut_slip=[centro(:,3)> -40000];
%  slipf=zeros(size(centro,1),1); % para cada triangulo
%  slipf(cut_slip)=1;
%%


%%
%gps=load('./plot_gmt/gps_inter.txt');
gps= load('gps_iqq.txt');


% %seleciona solo gps norte de -45S 
% cut=gps(:,2)<-45;
% gps(cut,:)=[];
% 
% cut=gps(:,1)>-62;
% gps(cut,:)=[];
% 
% save gps gps -ASCII


% gps en lon lat a coordenadas UTM
[x_gps,y_gps]=lat2utm_f(gps(:,1),gps(:,2)); %% en m
% coordenada vertical estaciones =0 m
z_gps=x_gps*0;

n_GPS=size(x_gps,1);
% 

%% % Create Green Function GF relacionan el moviento en la falla con mov. en GPS
%% 
%% hacerlo cuando cambie la distribucion de GPS 
%% Las GFs se recalculan solo si cambia distribucion de GPS o nueva falla
do_gf=2
if do_gf>1

h = waitbar(0,'Please wait...');

% Funciones de Green deslizamiento direccion dip
GX_ew=[]; GX_ns=[]; GX_up=[];

% Funciones de Green deslizamiento direccion strike
GY_ew=[];GY_ns=[]; GY_up=[];
 
%% calcula 
for i=1:size(tri,1)
%     i
  x_t=[xyz(tri(i,1),1) ; xyz(tri(i,2),1);  xyz(tri(i,3),1)];
  y_t=[xyz(tri(i,1),2) ; xyz(tri(i,2),2);  xyz(tri(i,3),2)];
  z_t=[xyz(tri(i,1),3) ; xyz(tri(i,2),3);  xyz(tri(i,3),3)];
 
% % guarda coordendas de los vertices para plot gmt en lat/lon
  [la1,lo1]=utm2lat_f(x_t(1),y_t(1));
  [la2,lo2]=utm2lat_f(x_t(2),y_t(2));
  [la3,lo3]=utm2lat_f(x_t(3),y_t(3));
%  
  vertices(i,:)=[lo1 la1 lo2 la2 lo3 la3];
  save vertices vertices -ASCII
  
 P1=[x_t(1) y_t(1) z_t(1)];  P2=[x_t(2) y_t(2) z_t(2)];
 P3=[x_t(3) y_t(3) z_t(3)];

 % % strike slip : deslizamiento lateral
  [ue_ss,un_ss,uv_ss] = TDdispHS(x_gps,y_gps,gps(:,2)*0,P1,P2,P3,1,0,0,.25);
% % dip
% -1 para cosismicos 1 para backslip
  [ue_ds,un_ds,uv_ds] = TDdispHS(x_gps,y_gps,gps(:,2)*0,P1,P2,P3,0,-1,0,.25);
  
% % dip slip deslizamiento vertical
   GX_ew=[GX_ew ue_ds];
   GX_ns=[GX_ns un_ds];
   GX_up=[GX_up uv_ds];
% % strike  
   GY_ew=[GY_ew ue_ss];
   GY_ns=[GY_ns un_ss];
   GY_up=[GY_up uv_ss];
% 
  waitbar(i/size(tri,1),h)

end

%une la matriz
G=[GX_ew GY_ew; GX_ns GY_ns; GX_up  GY_up];
save G.txt G -ASCII 
close(h)

else 
    G=load('G.txt');
    vertices=load('vertices');
end


%% estima desplazamiento total

% slip dip= slipf (slip_max) % slip_strike=0.25*slipf
%slip_dip= -slipf*1;         %cambia a negativa para q vay al oeste (falla inversa)
%slip_strike=0.25*slipf;   % 0.25 cm 

slip_dip= slipf* 8;       % terremoto -  %falla normal : + backslip=plate velocity
slip_strike= slipf* 0.5;

slip_total=sqrt(slip_dip.^2+slip_strike.^2);

save slip_dip.txt slip_dip -ASCII
save slip_strike.txt slip_strike -ASCII

%multiplicacion de G*Slip= desplazamiento
dtotal=G*[slip_dip;slip_strike];

%desplazamientos modelados
dew=dtotal(1:n_GPS);
dns=dtotal(1+n_GPS:n_GPS*2);
dup=dtotal(1+2*n_GPS:end);

% % plot displacements
figure(2); 
trisurf(tri,xyz(:,1),xyz(:,2),xyz(:,3),slip_total(:));
axis equal,colorbar
colormap jet
hold on
plot3(x_gps,y_gps,z_gps , 'o');

 plot3(x_gps,y_gps, y_gps *0,'or','MarkerSize',5)
 quiver3(x_gps,y_gps, y_gps *0,  dew ,dns , dup)
 ylim([min(y_gps), max(y_gps)])
   

%% estima la magnitud
 %% magnitude
    y=90e9; %young modul
    p=0.25 ; %possion
    s=0.5* (y /(1+p)); %shear modulus

    for i=1:size(slip_total,1)
    mo(i)=area(i)*slip_total(i); %area * slip  
    end

    moment=sum(mo)*s;
    magnitud=(2/3)*( log10(moment) - 9.1)
    
       
    fid=fopen('magnitude.txt','w');
    fprintf(fid,'%2.1f\n', magnitud);
    fclose(fid);
    
%%%%


%% cambia falla a long-lat para plot
[nodos_lat,nodos_lon]= utm2lat_f(coord(:,2),coord(:,3)); 
[centro_lat,centro_lon]=utm2lat_f(centro(:,1),centro(:,2));

%% plot para gmt
%% save triangulos para gmt
% para plot de dip slip
%save slip total

vertices_dip=[vertices ,slip_total];
save('slipx_tri.txt','vertices_dip','-ASCII');
system('./format_coords.py slipx_tri.txt');

out_slip_centro=[centro_lon,centro_lat,slip_total];
save out_slip.txt out_slip_centro -ASCII

% predicciones desplazamietno en GPS (errores =10 % de la
% prediccion) en mm
obs_h=[gps(:,1) , gps(:,2), dew, dns, dew*0.1,dns*0.1];
save obs.txt obs_h -ASCII
    
obs_v=[gps(:,1) , gps(:,2), dup*0, dup, dup*0.1 , dup*0.1 ];
save obs_up.txt obs_v -ASCII

 system('mv obs.txt out_slip.txt gmt_input.txt magnitude.txt  obs_up.txt ./plot_gmt');
% 
% cd plot_gmt
% system('./plot_slip.gmt')

%cd ..
