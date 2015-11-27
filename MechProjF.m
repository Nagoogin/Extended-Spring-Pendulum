function MechProjF(dt,T,l1s,r,l_2,thetai,k,m)
% Analytical Mechanics Group Project - Spring 2015
% Kevin Nguyen & Grant Young
% Double Pendulum (w/ Spring & single mass[m])

% Acceleration functions as the Lagrangian equations of motion for this
% particular system

% Call MechProj(.05,30,4,5,4,pi/4,8,1)

N       = round(T/dt); % setting arrays for all our desired values
t       = NaN(1,N);
theta   = NaN(1,N);
dtheta  = NaN(1,N);
dl1     = NaN(1,N);
l1      = NaN(1,N);
l2      = ones(1,N)*l_2;
l3      = NaN(1,N);
z       = zeros(1,N);
t2      = 0:.01:2*pi;

n = 1; % initialization 
t(n)      = 0;
theta(n)  = thetai;
dtheta(n) = 0;
l1(n)     = l1s;
dl1(n)    = 0;

while n<N
    dl1m = dl1(n)+ddl1(l1(n),r,l2(n),theta(n),dtheta(n),k,m)*0.5*dt;
    l1m  = l1(n)+dl1(n)*0.5*dt;
    
    dthetam = dtheta(n)+ddTheta(l1(n),dl1(n),l2(n),theta(n),dtheta(n))*0.5*dt;
    thetam  = theta(n)+dtheta(n)*0.5*dt;
    
    t(n+1)   = t(n)+dt;
    
    dl1(n+1) = dl1(n)+ddl1(l1m,r,l2(n),thetam,dthetam,k,m)*dt;
    l1(n+1)  = l1(n)+dl1m*dt;
    
    dtheta(n+1) = dtheta(n)+ddTheta(l1m,dl1m,l2(n),thetam,dthetam)*dt;
    theta(n+1)  = theta(n)+dthetam*dt;
    
    l3(n) = l1(n)+l2(n);
    
    figure(1) % plotting the motion within the while loop
    P = polar(t2,15*ones(size(t2))); % sets an invisible boundary (sets window of the polar plot)
    set(P,'Visible','off');
    view([-90,91]); % rotates the polar plot so the 0 degrees is straight down
    a = [z(n),theta(n)];
    b = [z(n),l1(n)];
    c = [theta(n),theta(n)];
    d = [l1(n),l3(n)];
    hold on
    polar(a,b,'g');
    polar(c,d,'r');
    polar(theta(n),l1(n));
    polar(theta(n),l3(n),'ok');
    polar(theta,l3,'c');
    hold off
    pause(0.0001);
    
    n = n+1;
end
end

% acceleration functions
function a1 = ddl1(l1,r,l2,theta,dtheta,k,m) 
g = 9.81;
a1 = (l1*(dtheta^2))+(l2*(dtheta^2))-((k/m)*(l1-r))+(g*cos(theta));
end

function a2 = ddTheta(l1,dl1,l2,theta,dtheta)
g = 9.81;
a2 = ((-g*sin(theta)*(l1+l2))/(2*l1*l2+l1^2+l2^2))-((2*dl1*dtheta*(l1+l2))/(2*l1*l2+l1^2+l2^2));
end