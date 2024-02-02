clc
clear

wrkDir = './' ;
% Loading data for dt = 0.16
dataStr = strcat(wrkDir,'x016.mat') ;
load(dataStr);

% Loading data for dt = 0.08
dataStr = strcat(wrkDir,'x008.mat') ;
load(dataStr);

% Loading data for dt = 0.04
dataStr = strcat(wrkDir,'x004.mat') ;
load(dataStr);

% Loading data for dt = 0.02
dataStr = strcat(wrkDir,'x002.mat') ;
load(dataStr);

% Loading data for dt = 0.01
dataStr = strcat(wrkDir,'x001.mat') ;
load(dataStr);

% Loading data for dt = 0.005
dataStr = strcat(wrkDir,'x0005.mat') ;
load(dataStr);


%calculating errors 
error1 = norm(T1 - T6)/norm(T6);
error2 = norm(T2 - T6)/norm(T6);
error3 = norm(T3 - T6)/norm(T6);
error4 = norm(T4 - T6)/norm(T6);
error5 = norm(T5 - T6)/norm(T6);

dt = [0.08,0.04,0.02,0.01];
error = [ error2, error3, error4, error5];

%plotting
figure(6)
loglog(dt,error);
hold on
scatter(dt,error);
xlabel("dt")
ylabel("Error")
title("Log Log graph of error vs dt")
legend("Slope is 2.1088")
xlim tight 
ylim tight

p = polyfit(log(dt),log(error),1) ;
slope = p(1);