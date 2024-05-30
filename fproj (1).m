clc;
clear all;
close all;
mu = 3.986e14;
RAAN = deg2rad(60);%deg
e = 0.2;
i =deg2rad(55);%deg
f=0;
AOP=deg2rad(53.2162);%deg
GMST =175.4942;%deg
T = 2*60*60; %sec
a=(T*sqrt(mu)/(2*pi))^(2/3);

%% 1
r_eci = [-634.621106e03;4824.19897e03;4229.74236e03];  %analytical solution
v_eci = [-6.01121473e03;-4.49508378e03;4.22492277e03]; %analytical solution
x0 = [r_eci,v_eci];
tspan = linspace(0,T,1000);
options = odeset('RelTol',1e-6,'AbsTol',1e-4);
[t,x] = ode45(@(t,x)TwoBP(t,x,mu),tspan,x0,options);

xp = x(:,1);
yp = x(:,2);
zp = x(:,3);
figure
plot3(xp,yp,zp);
title('ECI Frame');
xlabel ('x');
ylabel ('y');
zlabel ('z');
grid on;

%% 2
%ECI to ECEF
for i=1:length(x)
    [r_ecef(i,:), v_ecef(i,:)] = ECI2ECEF(x(i,1:3),x(i,4:6),AOP,T,GMST,mu); 
end

xps = r_ecef(:,1);
yps = r_ecef(:,2);
zps = r_ecef(:,3);
figure
plot3(xps,yps,zps);
title('ECEF Frame');
xlabel ('x');
ylabel ('y');
zlabel ('z');
grid on;

figure 
plot(xps,yps);
title('X-Y plane for ECEF Frame');
xlabel('X axis');
ylabel('Y axis');
figure 
plot(yps,zps);
title('Y-Z plane for ECEF Frame');
xlabel('Y axis');
ylabel('Z axis');
figure 
plot(xps,zps);
title('X-Z plane for ECEF Frame');
xlabel('X axis');
ylabel('Z axis');




%% 3
% Convert ECEF to latitude, longitude
for i =1:length(xps)
    lon(i)= atan2d(yps(i),xps(i));
    lat(i) = atan2d(zps(i),sqrt(xps(i)^2+yps(i)^2));
end

figure
load('topo.mat','topo');
topoplot = [topo(:,181:360),topo(:,1:180)];
contour(-180:179,-90:89,topoplot,[0,0],'black','linewidth',1);
grid on
grid minor
axis equal
hold on
plot(lon',lat','m.','LineWidth',2)
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
set(gca,'FontSize',18)

%% 4
rp = [-634.621106e03; 4824.19897e03; 4229.74236e03];
vp = [-6.01121473e03; -4.49508378e03; 4.22492277e03];
T = 7200;
tspan = linspace(0, T, 1000);
for i = 2:length(tspan)
    [r1(:, i), v1(:, i)] = FGfunc(rp, vp, tspan(i), mu);
end
figure
plot3(r1(1,:), r1(2,:), r1(3,:),'m.');
title('F & G Propogation')
xlabel ('x');
ylabel ('y');
zlabel ('z');


%% 5 Position
figure
plot3(xp,yp,zp);
hold on;
plot3(r1(1,:), r1(2,:), r1(3,:),'m.');
hold on;
plot3(xps,yps,zps);
title('Comparison');
legend('ECI','F&G','ECEF');
xlabel ('x');
ylabel ('y');
zlabel ('z');
%% 5 Velocity
figure
plot3(x(:,4),x(:,5),x(:,6));
hold on;
plot3(v1(1,:), v1(2,:), v1(3,:));
hold on;
plot3(v_ecef(:,1), v_ecef(:,2), v_ecef(:,3));
title('Comparison');
legend('ECI','F&G','ECEF');
xlabel ('x');
ylabel ('y');
zlabel ('z');

%% 6  in ECI Frame
rf = [x(1000,1);x(1000,2);x(1000,3)]; % ECI
vf = [x(1000,4);x(1000,5);x(1000,6)]; %ECI 

[a_final, e_final, i_final, RAAN_final, AOP_final, f_final] = RV2OE(rf,vf,mu)


%%
%% Functions
%%

%%  Function 1: FGFunc
function [r,v] = FGfunc(r0,v0,dt,mu)
delta = 1e-2;
    [anew, e1, ~, ~, ~, f0] = RV2OE(r0, v0, mu);
    E0=atan(sqrt((1-e1)/(1+e1))*tand(f0/2))*2;
    M0=E0-(e1*sin(E0));
    E = kepler(anew, e1, mu, dt,0, delta);
    delaE = E - M0 ;
    [r, v] = fg_eqs(anew, delaE, r0, v0, norm(r0), mu,dt);

function E = kepler(a,e,mu,t,tp,delta)
        M=(sqrt(mu/a^3)*(t-tp));
        eo=M;
        en=0;
        error=2*delta;
        while(error>=delta)
            en=eo-(((M+e*sin(eo))-eo)/((e*cos(eo))-1));
            error=abs(en-eo);
            eo=en;
        end
        E=en;
    end
    function [r, v] = fg_eqs(a, delaE, r0, v0, rma, mu,dt)
        f = 1 - (a/rma * (1 - cos(delaE)));
        g =  dt-( sqrt(a^3/mu) * (delaE - sin(delaE)));
        r = (f * r0) + (g * v0);
        rc= norm(r);
        fdot = -sqrt(mu*a/(rma*rc)) * sin(delaE);
        gdot = 1 - (a/rc * (1 - cos(delaE)));
        v = fdot * r0 + gdot * v0;
    end
end
%% Function 2: TwoBP
function dx = TwoBP(t, x, mu)
    r = x(1:3);
    v = x(4:6);
    a = -mu * r / norm(r)^3;
    dx = [v;a];
end
%% Function 3: RV2OE
function [a, e, i, RAAN, AOP, f] = RV2OE(r, v, mu)
    r_mag = norm(r);
    v_mag = norm(v);
    h = cross(r, v);
    h_mag = norm(h);
    N = cross([0; 0 ;1], h);
    N_mag = norm(N);
    e_vec = (1/mu) * ((v_mag^2 - mu/r_mag).*r - dot(r, v).*v);
    e = norm(e_vec);
    a = 1 ./ (2./r_mag - v_mag.^2./mu);
    i = acos(h(3) / h_mag);
    RAAN = acos(N(1) / N_mag);
    if N(2) < 0
        RAAN = 2*pi - RAAN;
    end
    AOP = acos(dot(N, e_vec) / (N_mag * e));
    if e_vec(3) < 0
        AOP = 2*pi - AOP;
    end
    f = acos(dot(e_vec, r) / (e * r_mag));
    if dot(r, v) < 0
        f = 2*pi - f;
    end
    i = rad2deg(i);
    RAAN = rad2deg(RAAN);
    AOP = rad2deg(AOP);
    f = rad2deg(f);
end


%% Function 4: OE2RV
function [r,v] = OE2RV(a,e,I,RAAN,AOP,f,mu)

    p = a*(1 - e^2);
    r = (p/(1 + e*cos(f)))*[cos(f); sin(f); 0];
    v = sqrt(mu/p)*[-sin(f); e + cos(f); 0];

    theta = AOP+f;
    
    a = [(cos(RAAN)*cos(theta))-(sin(RAAN)*sin(theta)*sin(I)),(sin(RAAN)*cos(theta))+(cos(RAAN)*sin(theta)*cos(I)), sin(RAAN)*sin(I);...
        -(cos(RAAN)*sin(theta))-(sin(RAAN)*cos(theta)*cos(I)),-(sin(RAAN)*cos(theta))+(cos(RAAN)*cos(theta)*cos(I)), cos(RAAN)*sin(I);...
    sin(RAAN)*sin(I), -cos(RAAN)*sin(I), cos(I)];
    r= a'*r;
    v=a'*v;
end
%% Function 5: ECI2ECEF
function [r_ecef,v_ecef] = ECI2ECEF(r_eci,v_eci,omega,t,GMST,mu)

    theta = GMST + omega*t;

    r3_theta = [cos(theta), -sin(theta), 0; 
        sin(theta), cos(theta), 0;
        0, 0, 1];

    r_ecef = r3_theta*r_eci';
    v_ecef = r3_theta*v_eci';
end
