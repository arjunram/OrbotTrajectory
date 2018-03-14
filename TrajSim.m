clear 
clc
k = 1;
tf = 7;
tspan = [0 tf];
Xc0 = [10; 2; 90];
syms T;
xr = 5*cos(T);
yr = 5*sin(T);
fplot(xr,yr,[0 tf])
xrdot = diff(xr);
yrdot = diff(yr);
theta = atan2(yrdot,xrdot);
Xr(1,1) = xr;
Xr(2,1) = yr;
Xr(3,1) = theta;
Xrdot(1,1) = xrdot;
Xrdot(2,1) = yrdot;
Xrdot(3,1) = (xrdot*diff(yrdot) - yrdot*diff(xrdot))/(xrdot^2 + yrdot^2);
[t,Xc] = ode45(@(t,Xc) odefun(t,Xr,Xrdot,Xc,k),tspan,Xc0);
hold on;
plot(Xc(:,1),Xc(:,2),'r');

Xcsim = Xc0;
scatter(Xcsim(1),Xcsim(2));
Tnew = 0:tf/500:tf;
for i = 2:length(Tnew)
    Xcdot = odefun(Tnew(i),Xr,Xrdot,Xcsim,k);
    Xcsim = Xcsim + Xcdot*(Tnew(i)-Tnew(i-1));
    scatter(Xcsim(1),Xcsim(2),'b');
end
xlim([-6 6])
ylim([-6 6])
hold off

function Xcdot = odefun(t,Xr,Xrdot,Xc,k)
c = cos(Xc(3));
s = sin(Xc(3));
R1 = [c -s 0; s c 0;0 0 1];
R2 = [c s -c*(Xr(2)-Xc(2))+s*(Xr(1)-Xc(1)); -s c c*(Xr(1)-Xc(1))+s*(Xr(2)-Xc(2));0 0 1];
z = inv(R1)*inv(R2)*inv(R1)*Xrdot + k*inv(R1)*inv(R2)*(Xr-Xc);
Xcdot = double(subs(Xrdot - inv(R1)*Xrdot + R2*R1*z,t));
end

function z = inputs(t,Xr,Xrdot,Xc,k)
c = cos(Xc(3));
s = sin(Xc(3));
R1 = [c -s 0; s c 0;0 0 1];
R2 = [c s -c*(Xr(2)-Xc(2))+s*(Xr(1)-Xc(1)); -s c c*(Xr(1)-Xc(1))+s*(Xr(2)-Xc(2));0 0 1];
z = double(subs(inv(R1)*inv(R2)*inv(R1)*Xrdot + k*inv(R1)*inv(R2)*(Xr-Xc),t));
%Xcdot = double(subs(Xrdot - inv(R1)*Xrdot + R2*R1*z,t));
end
