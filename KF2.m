clear all,clc
% Basic kalman filter simulation
% Author: Yunxi TANG
% Date: 21/08/2019

% state space transfer function:
% s(k)=s(k-1)+v(k-1)*dt+0.5*a*dt*dt
% v(k)=v(k-1)+a*dt;
% matrix form:
% x(k) = A*x(k-1)+B*u(k)+w(k)
% x(k) = [s(k);
%        v(k)]
% A = [1 dt;
%      0 1]
% B = [0.5*dt*dt;
%     dt]
% w(k)~N(0,Q) --> noise

% Measure Method in observe space
% z(k) = H(k)x(k)+v(k)
% z(k) --> measured vector
% H(k) --> observe matrix   
% v(k)~N(0,R) --> noise

m = 5; F = 2;
a = F / m;

dt = 0.1;
t = 0:dt:10;
len = length(t);

A = [1 dt;
     0 1];
B = [0.5*dt*dt;
     dt];
 
H = [1 0;
     0 1];
 
q = 1;
r = 1.5;

% experimental data
w = 0.0 + q.*randn(2,len+1,'double');
v = 0.0 + r.*randn(2,len+1,'double');

% covriance matrix
Q = [q*q 0;
     0 q*q];
R = [r^2 0;
     0 r^2];

% initialize P and K
P = [0 0;
     0 0];
K = [0 0;
     0 0];

%% initial state                 
z(:,1) = [0;0];                             %guancezhi

%% ideal condition x_ide, prediction x and sensor data z
x_ide(:,1) = [0;0];
z(:,1) = [0;0];
x(:,1) = [0;0];
for j = 2:len
    x_ide(:,j) = A*x_ide(:,j-1)+B*a;        %ide
    x(:,j) = A*x(:,j-1) + w(:,j-1);         %pre
    z(:,j) = H*x(:,j)+v(:,j);               % sensor data
end

%% filter
x_pre(:,1) = [0;0];

for k = 2:len
    x_pre(:,k) = A*x_pre(:,k-1) + B*a + w(:,k-1);    % one step prediction
    
    P = A*P*A'+Q;
    
    K = P*H'\(H*P*H'+R);
    
    x_pre(:,k) = x_pre(:,k) + K*(z(:,k)-H*x_pre(:,k-1));
    
    P = P - K*H*P;  
    
end

%% Figure
figure(1);
plot(t,x_pre(1,1:len),'-.',t,z(1,1:len),t,x_ide(1,1:len),'--');
legend('Filter','Observation','Prediction');
title('s');

figure(2);
plot(t,x_pre(2,1:len),'-.',t,z(2,1:len),t,x_ide(2,1:len),'--');
legend('Filter','Observation','Prediction');
title('V');


    
    

