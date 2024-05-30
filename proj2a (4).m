clear all
clc
%% Part B- Testing the alogorithm

a_semi1=13000;
a_semi2=7226.58;
mu = 3.986e5;
e1=0.3;
e2=0.444819;
p1=a_semi1*(1-(e1^2));
p2=a_semi2*(1-(e2^2));
aop1=deg2rad(50);
aop2=deg2rad (301.901);
RAAN1=deg2rad(30);
I1=deg2rad(20);
h1=sqrt(mu*a_semi1*(1-e1^2));
h2=sqrt(mu*a_semi2*(1-e2^2));
deltaaop=aop2-aop1;
alpha = (e2*cos(deltaaop))-e1;
beta =(e1*p2)-(e2*p1*cos(deltaaop));
gamma = e1*e2*sin(deltaaop);

a=((e1^2-1)/(e1^2))-(alpha^2/gamma^2);
b=((2*p1)/e1^2)-((2*alpha*beta)/gamma^2);
c=-((p1^2/e1^2)+(beta^2/gamma^2));

r1= (-b+sqrt(b^2-(4*a*c)))/(2*a);
r2= (-b-sqrt(b^2-(4*a*c)))/(2*a);
f1=acos(((p1/r1)-1)/e1);
f2= acos(((p2/r2)-1)/e2);
theta1=aop1;
theta2=aop2;

r1vec=intialcon1(r1,f1);
a1dcm=dcm(RAAN1,theta1,I1);
x1=a1dcm*r1vec;
r2vec=intialcon1(r2,f2);
a2dcm=dcm(RAAN1,theta2,I1);
x2=a2dcm*r2vec;
v1vec= intialcon2(h1,f1,mu,e1);
vx1=a1dcm*v1vec;
v2vec= intialcon2(h2,f2,mu,e2);
vx2=a2dcm*v2vec;
deltav1=abs(norm(vx2)-norm(vx1))
%% Numerical Part
disp('Numerical')
%%
load('IODMeasurements2.mat');
mu = 3.986e+05;
%%
disp('Orbit 1')
[r1,r2,r3] =gauss(1);
[r_1,v_1]= gibbs(r1,r2,r3)
[a_1, e_1, i_1, omega_1, w_1, f_1]=RV2OE(r_1,v_1,mu)
T =2*pi*sqrt(a_1^3/mu);
x0 = [r_1,v_1];
tspan = linspace(0,T,1000);
options = odeset('RelTol',1e-6,'AbsTol',1e-4);
[t,x] = ode45(@(t,x)TwoBP(t,x,mu),tspan,x0,options);
xp = x(:,1);
yp = x(:,2);
zp = x(:,3);
 
disp('Orbit 2')

[r1,r2,r3] =gauss(4);
[r_2,v_2]= gibbs(r1,r2,r3);
[a_2, e_2, i_2, omega_2, w_2, f_2]=RV2OE(r_2,v_2,mu)
T =2*pi*sqrt(a_1^3/mu);
x0 = [r_2,v_2];
tspan = linspace(0,T,1000);
options = odeset('RelTol',1e-6,'AbsTol',1e-4);
[t,x] = ode45(@(t,x)TwoBP(t,x,mu),tspan,x0,options);
xp2 = x(:,1);
yp2 = x(:,2);
zp2 = x(:,3);


r1vec=intialcon1(norm(r_1),f_1);
a1dcm=dcm(omega_1,w_1,i_1);
x1=a1dcm*r1vec;
r2vec=intialcon1(norm(r_2),f_2);
a2dcm=dcm(omega_2,w_2,i_2);
x2=a2dcm*r2vec;
h1=sqrt(mu*a_1*(1-e_1^2));
h2=sqrt(mu*a_2*(1-e_2^2));
v1vec= intialcon2(h1,f_1,mu,e_1);
vx1=a1dcm*v1vec;
v2vec= intialcon2(h2,f_2,mu,e_2);
vx2=a2dcm*v2vec;
deltav=abs(norm(vx2)-norm(vx1))

%%
disp('Orbit 3')
[r1,r2,r3] =gauss(7);
[r_1,v_1]= gibbs(r1,r2,r3);
[a_1, e_1, i_1, omega_1, w_1, f_1]=RV2OE(r_1,v_1,mu)
T =2*pi*sqrt(a_1^3/mu);
x0 = [r_1,v_1];
tspan = linspace(0,T,1000);
options = odeset('RelTol',1e-6,'AbsTol',1e-4);
[t,x] = ode45(@(t,x)TwoBP(t,x,mu),tspan,x0,options);
xp3 = x(:,1);
yp3 = x(:,2);
zp3 = x(:,3);
r1vec=intialcon1(norm(r_1),f_1);
a1dcm=dcm(omega_1,w_1,i_1);
x1=a1dcm*r1vec;
r2vec=intialcon1(norm(r_2),f_2);
a2dcm=dcm(omega_2,w_2,i_2);
x2=a2dcm*r2vec;
h1=sqrt(mu*a_1*(1-e_1^2));
h2=sqrt(mu*a_2*(1-e_2^2));
v1vec= intialcon2(h1,f_1,mu,e_1);
vx1=a1dcm*v1vec;
v2vec= intialcon2(h2,f_2,mu,e_2);
vx2=a2dcm*v2vec;
deltav=abs(norm(vx2)-norm(vx1))

disp('Orbit 4')
[r1,r2,r3] =gauss(10);
[r_2,v_2]= gibbs(r1,r2,r3);
[a_2, e_2, i_2, omega_2, w_2, f_2]=RV2OE(r_2,v_2,mu)
T =2*pi*sqrt(a_1^3/mu);
x0 = [r_2,v_2];
tspan = linspace(0,T,1000);
options = odeset('RelTol',1e-6,'AbsTol',1e-4);
[t,x] = ode45(@(t,x)TwoBP(t,x,mu),tspan,x0,options);
xp4 = x(:,1);
yp4 = x(:,2);
zp4 = x(:,3);

r1vec=intialcon1(norm(r_1),f_1);
a1dcm=dcm(omega_1,w_1,i_1);
x1=a1dcm*r1vec;
r2vec=intialcon1(norm(r_2),f_2);
a2dcm=dcm(omega_2,w_2,i_2);
x2=a2dcm*r2vec;
h1=sqrt(mu*a_1*(1-e_1^2));
h2=sqrt(mu*a_2*(1-e_2^2));
v1vec= intialcon2(h1,f_1,mu,e_1);
vx1=a1dcm*v1vec;
v2vec= intialcon2(h2,f_2,mu,e_2);
vx2=a2dcm*v2vec;
deltav=abs(norm(vx2)-norm(vx1))

%%
disp('Orbit 5')
[r1,r2,r3] =gauss(13);
[r_1,v_1]= gibbs(r1,r2,r3);
[a_1, e_1, i_1, omega_1, w_1, f_1]=RV2OE(r_1,v_1,mu)
T =2*pi*sqrt(a_1^3/mu);
x0 = [r_1,v_1];
tspan = linspace(0,T,1000);
options = odeset('RelTol',1e-6,'AbsTol',1e-4);
[t,x] = ode45(@(t,x)TwoBP(t,x,mu),tspan,x0,options);
xp5 = x(:,1);
yp5 = x(:,2);
zp5 = x(:,3);
r1vec=intialcon1(norm(r_1),f_1);
a1dcm=dcm(omega_1,w_1,i_1);
x1=a1dcm*r1vec;
r2vec=intialcon1(norm(r_2),f_2);
a2dcm=dcm(omega_2,w_2,i_2);
x2=a2dcm*r2vec;
h1=sqrt(mu*a_1*(1-e_1^2));
h2=sqrt(mu*a_2*(1-e_2^2));
v1vec= intialcon2(h1,f_1,mu,e_1);
vx1=a1dcm*v1vec;
v2vec= intialcon2(h2,f_2,mu,e_2);
vx2=a2dcm*v2vec;
deltav=abs(norm(vx2)-norm(vx1))


disp('Orbit 6')
[r1,r2,r3] =gauss(16);
[r_2,v_2]= gibbs(r1,r2,r3);
[a_2, e_2, i_2, omega_2, w_2, f_2]=RV2OE(r_2,v_2,mu)
T =2*pi*sqrt(a_2^3/mu);
x0 = [r_2,v_2];
tspan = linspace(0,T,1000);
options = odeset('RelTol',1e-6,'AbsTol',1e-4);
[t,x] = ode45(@(t,x)TwoBP(t,x,mu),tspan,x0,options);
xp6 = x(:,1);
yp6 = x(:,2);
zp6 = x(:,3);

r1vec=intialcon1(norm(r_1),f_1);
a1dcm=dcm(omega_1,w_1,i_1);
x1=a1dcm*r1vec;
r2vec=intialcon1(norm(r_2),f_2);
a2dcm=dcm(omega_2,w_2,i_2);
x2=a2dcm*r2vec;
h1=sqrt(mu*a_1*(1-e_1^2));
h2=sqrt(mu*a_2*(1-e_2^2));
v1vec= intialcon2(h1,f_1,mu,e_1);
vx1=a1dcm*v1vec;
v2vec= intialcon2(h2,f_2,mu,e_2);
vx2=a2dcm*v2vec;
deltav=abs(norm(vx2)-norm(vx1))



%%
disp('Orbit 7')
[r1,r2,r3] =gauss(19);
[r_1,v_1]= gibbs(r1,r2,r3);
[a_1, e_1, i_1, omega_1, w_1, f_1]=RV2OE(r_1,v_1,mu)
T =2*pi*sqrt(-a_1^3/mu);
x0 = [r_1,v_1];
tspan = linspace(0,T,1000);
options = odeset('RelTol',1e-6,'AbsTol',1e-4);
[t,x] = ode45(@(t,x)TwoBP(t,x,mu),tspan,x0,options);
xp6 = x(:,1);
yp6 = x(:,2);
zp6 = x(:,3);


r1vec=intialcon1(norm(r_1),f_1);
a1dcm=dcm(omega_1,w_1,i_1);
x1=a1dcm*r1vec;
r2vec=intialcon1(norm(r_2),f_2);
a2dcm=dcm(omega_2,w_2,i_2);
x2=a2dcm*r2vec;
h1=sqrt(mu*a_1*(1-e_1^2));
h2=sqrt(mu*a_2*(1-e_2^2));
v1vec= intialcon2(h1,f_1,mu,e_1);
vx1=a1dcm*v1vec;
v2vec= intialcon2(h2,f_2,mu,e_2);
vx2=a2dcm*v2vec;
deltav=abs(norm(vx2)-norm(vx1))

disp('Orbit 8')
[r1,r2,r3] =gauss(22);
[r_2,v_2]= gibbs(r1,r2,r3);
[a_2, e_2, i_2, omega_2, w_2, f_2]=RV2OE(r_2,v_2,mu)
T =2*pi*sqrt(-a_2^3/mu);
x0 = [r_2,v_2];
tspan = linspace(0,T,1000);
options = odeset('RelTol',1e-6,'AbsTol',1e-4);
[t,x] = ode45(@(t,x)TwoBP(t,x,mu),tspan,x0,options);
xp7 = x(:,1);
yp7 = x(:,2);
zp7 = x(:,3);
r1vec=intialcon1(norm(r_1),f_1);
a1dcm=dcm(omega_1,w_1,i_1);
x1=a1dcm*r1vec;
r2vec=intialcon1(norm(r_2),f_2);
a2dcm=dcm(omega_2,w_2,i_2);
x2=a2dcm*r2vec;
h1=sqrt(mu*a_1*(1-e_1^2));
h2=sqrt(mu*a_2*(1-e_2^2));
v1vec= intialcon2(h1,f_1,mu,e_1);
vx1=a1dcm*v1vec;
v2vec= intialcon2(h2,f_2,mu,e_2);
vx2=a2dcm*v2vec;
deltav=abs(norm(vx2)-norm(vx1))

figure
plot3(xp,yp,zp);
hold on
plot3(xp3,yp3,zp3);
plot3(xp2,yp2,zp2);
plot3(xp4,yp4,zp4);
plot3(xp5,yp5,zp5);
%plot3(xp6,yp6,zp6);
%plot3(xp7,yp7,zp7);


xlabel ('x');
ylabel ('y');
zlabel ('z');
grid on;

% [x1,y1,z1]= plot_orbit(a_1,e_1,i_1,omega_1,w_1);
% [r1,r2,r3] =gauss(4);
% [r_2,v_2]= gibbs(r1,r2,r3,4);
% [a_2, e_2, i_2, omega_2, w_2, f_2]=RV2OE(r_2,v_2,mu);
% [x2,y2,z2]= plot_orbit(a_2,e_2,i_2,omega_2,w_2);
% [r1,r2,r3] =gauss(7);
% [r_3,v_3]= gibbs(r1,r2,r3,7);
% [a_3, e_3, i_3, omega_3, w_3, f_3]=RV2OE(r_3,v_3,mu);
% [x3,y3,z3]= plot_orbit(a_3,e_3,i_3,omega_3,w_3);
% [r1,r2,r3] =gauss(10);
% [r_4,v_4]= gibbs(r1,r2,r3,10);
% [a_4, e_4, i_4, omega_4, w_4, f_4]=RV2OE(r_4,v_4,mu);
% [x4,y4,z4]= plot_orbit(a_4,e_4,i_4,omega_4,w_4);
% [r1,r2,r3] =gauss(13);
% [r_5,v_5]= gibbs(r1,r2,r3,13);
% [a_5, e_5, i_5, omega_5, w_5, f_5]=RV2OE(r_5,v_5,mu);
% [x5,y5,z5]= plot_orbit(a_5,e_5,i_5,omega_5,w_5);
% [r1,r2,r3] =gauss(16);
% [r_6,v_6]= gibbs(r1,r2,r3,16);
% [a_6, e_6, i_6, omega_6, w_6, f_6]=RV2OE(r_6,v_6,mu);
% [x6,y6,z6]= plot_orbit(a_6,e_6,i_6,omega_6,w_6);
% [r1,r2,r3] =gauss(19);
% [r_7,v_7]= gibbs(r1,r2,r3,19);
% [a_7, e_7, i_7, omega_7, w_7, f_7]=RV2OE(r_7,v_7,mu);
% [x7,y7,z7]= plot_orbit(a_7,e_7,i_7,omega_7,w_7);
% [r1,r2,r3] =gauss(22);
% [r_8,v_8]= gibbs(r1,r2,r3,22);
% [a_8, e_8, i_8, omega_8, w_8, f_8]=RV2OE(r_8,v_8,mu);
% [x8,y8,z8]= plot_orbit(a_8,e_8,i_8,omega_8,w_8);
% 
%     figure;
%     plot3(x1, y1, z1);
%     hold on;
%     plot3(x2, y2, z2);
%     hold on;
%     plot3(x3, y3, z3);
%     hold on;
%     plot3(x4, y4, z4);
%     hold on;
%     plot3(x5, y5, z5);
%     hold on;
%     plot3(x6, y6, z6);
%     hold on;
%     plot3(x7, y7, z7);
%     hold on;
%     plot3(x8, y8, z8);
%     hold on;
%     xlabel('X Coordinate');
%     ylabel('Y Coordinate');
%     zlabel('Z Coordinate');
%     title('3D Orbital Plot');
%     grid on;
% % Define the array of inclination values
% inclinations = [1, 4, 7, 10, 13, 16];
% 
% % Define the standard gravitational parameter (assumed for Earth)
% mu = 3.986e+05;
% 
% % Prepare to plot all orbits in a single 3D plot
% figure;
% hold on;
% 
% % Loop over each inclination value
% for idx = 1:length(inclinations)
%     % Calculate position vectors using the Gauss method for current inclination
%     [r1, r2, r3] = gauss(inclinations(idx));
% 
%     % Calculate velocity vectors using the Gibbs method
%     [r, v] = gibbs(r1, r2, r3);
% 
%     % Convert position and velocity vectors to orbital elements
%     [a, e, i, omega, w, f] = RV2OE(r, v, mu);
% 
%     % Plot the orbit based on calculated orbital elements
%     [x, y, z] = plot_orbit(a, e, i, omega, w);
% 
%     % Add the current orbit plot to the 3D figure
%     plot3(x, y, z);
% end
% 
% % Customize the final plot
% xlabel('X Coordinate');
% ylabel('Y Coordinate');
% zlabel('Z Coordinate');
% title('3D Orbital Plot');
% grid on;
% hold off;

%% Function Gauss Method
function [r_1,r_2,r_3] = gauss(i)
load('IODMeasurements2.mat');
mu = 3.986e+05;

t1 = TIMES(i)-TIMES(i+1);
t3 = TIMES(i+2)-TIMES(i+1);
t13 = t3-t1;
%t13 = TIMES(i+2)-TIMES(i);

l1= [cosd(AZIMUTH(i))*cosd(ELEVATION(i)); 
    sind(AZIMUTH(i))*cosd(ELEVATION(i));
    sind(ELEVATION(i))];
l2=[cosd(AZIMUTH(i+1))*cosd(ELEVATION(i+1)); 
    sind(AZIMUTH(i+1))*cosd(ELEVATION(i+1));
    sind(ELEVATION(i+1))];
l3=[cosd(AZIMUTH(i+2))*cosd(ELEVATION(i+2)); 
    sind(AZIMUTH(i+2))*cosd(ELEVATION(i+2));
    sind(ELEVATION(i+2))];
r1= RSITES(i,:);
r2= RSITES(i+1,:);
r3= RSITES(i+2,:);

d0 = dot(l1,cross(l2,l3));

d11 = dot(r1', cross(l2,l3));
d12 = dot(r1', cross(l1,l3));
d13 = dot(r1', cross(l1,l2));

d21 = dot(r2', cross(l2,l3));
d22 = dot(r2', cross(l1,l3));
d23 = dot(r2', cross(l1,l2));

d31 = dot(r3', cross(l2,l3));
d32 = dot(r3', cross(l1,l3));
d33 = dot(r3', cross(l1,l2));

A = (1/d0) * (((-t3/t13)*d12)+d22+((t1/t13)*d32));
B = (1/(6*d0)) * ((-(t13^2-t3^2)*(t3/t13)*d12)+((t13^2-t1^2)*(t1/t13)*d32));

a = -(A^2) - (2*A*dot(l2,r2))-((norm(r2))^2);
b = -2*mu*B*(A+dot(l2,r2));
c= (-(mu)^2)*(B)^2;

eq = [1,0,a,0,0,b,0,0,c];

r=roots(eq);
r = r(r == real(r));
if r(1) > 0
   r = r(1);
elseif r(2) > 0
    r = r(2);
end
if(length(r)==1)
rho1 = (1/d0)*((6*((d31*(t1/t3))+(d21*(t13/t3)))*r^3 + mu*d31*(t13^2 - t1^2)*(t1/t3))/(6*r^3 + mu*(t13^2 - t3^2)) - d11);
rho2 = A + mu*B*(r^-3) ;
rho3 = (1/d0)*((6*((d13*(t3/t1))+(d23*(t13/t1)))*r^3 + mu*d13*(t13^2 - t3^2)*(t3/t1))/(6*r^3 + mu*(t13^2 - t1^2)) - d33);

r_1= r1' + (rho1*l1);
r_2= r2' + (rho2*l2);
r_3= r3' + (rho3*l3);
else 
rho1 = (1/d0)*((6*((d31*(t1/t3))+(d21*(t13/t3)))*r(2)^3 + mu*d31*(t13^2 - t1^2)*(t1/t3))/(6*r(2)^3 + mu*(t13^2 - t3^2)) - d11);
rho2 = A + mu*B*(r(2)^-3) ;
rho3 = (1/d0)*((6*((d13*(t3/t1))+(d23*(t13/t1)))*r(2)^3 + mu*d13*(t13^2 - t3^2)*(t3/t1))/(6*r(2)^3 + mu*(t13^2 - t1^2)) - d33);

r_1= r1' + (rho1*l1);
r_2= r2' + (rho2*l2);
r_3= r3' + (rho3*l3);

end
end

%% Gibbs Method 

function [pos,velo] = gibbs(r1,r2,r3)
mu = 3.986e+05;


eps= dot((r1/norm(r1)),(cross(r2,r3)/norm(cross(r2,r3))));

n = (norm(r1).*cross(r2,r3))+(norm(r2).*cross(r3,r1))+(norm(r3).*cross(r1,r2));
d = cross(r1,r2)+cross(r2,r3)+cross(r3,r1);
s = (r1.*(norm(r2)- norm(r3)))+(r2.*(norm(r3)- norm(r1)))+(r3.*(norm(r1)- norm(r2)));


velo = sqrt(mu/(dot(norm(n),norm(d))))*(((cross(d,r2)/norm(r2)))+s);
pos = r2;

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
%% Function 2: TwoBP
function dx = TwoBP(t, x, mu)
    r = x(1:3);
    v = x(4:6);
    a = -mu * r / norm(r)^3;
    dx = [v;a];
end
%% function delta v

function rp = intialcon1(rmag,fanal)
    rp =[rmag*cos(fanal);rmag*sin(fanal);0];
end
function vp = intialcon2(hfun,fanal1,mu1,efun)
 mubyh1=mu1/hfun;
 vp=mubyh1*[-sin(fanal1);(efun+cos(fanal1));0];
end
function adcm = dcm(RAAN,theta,I)
adcm = [(cos(RAAN)*cos(theta))-(sin(RAAN)*sin(theta)*cos(I)),(-sin(theta)*cos(RAAN))-(sin(RAAN)*cos(theta)*cos(I)), sin(RAAN)*sin(I);...
        (sin(RAAN)*cos(theta))+(cos(RAAN)*sin(theta)*cos(I)),-(sin(RAAN)*sin(theta))+(cos(RAAN)*cos(theta)*cos(I)), -cos(RAAN)*sin(I);...
    sin(theta)*sin(I),  cos(theta)*sin(I), cos(I)];
end