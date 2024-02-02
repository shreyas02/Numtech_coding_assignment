clc
clear

wrkDir = './' ;
% Loading data for dx = 0.1
dataStr = strcat(wrkDir,'x01.mat') ;
load(dataStr);

% Loading data for dx = 0.05
dataStr = strcat(wrkDir,'x005.mat') ;
load(dataStr);

% Loading data for dx = 0.025
dataStr = strcat(wrkDir,'x0025.mat') ;
load(dataStr);

% Loading data for dx = 0.0125
dataStr = strcat(wrkDir,'x00125.mat') ;
load(dataStr);

% Loading data for dx = 0.00625
dataStr = strcat(wrkDir,'x000625.mat') ;
load(dataStr);

%coordinates array for all mesh sizes 
elem1 = 1/0.1;
node1 = elem1 + 1;
crd1 = 1:node1;
crd1 = (crd1-1)/elem1;

elem2 = 1/0.05;
node2 = elem2 + 1;
crd2 = 1:node2;
crd2 = (crd2-1)/elem2;

elem3 = 1/0.025;
node3 = elem3 + 1;
crd3 = 1:node3;
crd3 = (crd3-1)/elem3;

elem4 = 1/0.0125;
node4 = elem4 + 1;
crd4 = 1:node4;
crd4 = (crd4-1)/elem4;

elem5 = 1/0.00625;
node5 = elem5 + 1;
crd5 = 1:node5;
crd5 = (crd5-1)/elem5;

%mapping each temp matrix 
for i=1:node1
    for j = 1:node1
        T11(i,j) = T1(i,j);
    end
end

for i=1:node2
    for j = 1:node2
        T22(i,j) = T2(i,j);
    end
end

for i=1:node3
    for j = 1:node3
        T33(i,j) = T3(i,j);
    end
end

for i=1:node4
    for j = 1:node4
        T44(i,j) = T4(i,j);
    end
end

for i=1:node5
    for j = 1:node5
        T55(i,j) = T5(i,j);
    end
end

for i=1:node1
    for j = 1:node1
        T51(i,j) = T5(16*(i-1) + 1,16*(j-1) + 1);
    end
end

for i=1:node2
    for j = 1:node2
        T52(i,j) = T5(8*(i-1) + 1,8*(j-1) + 1);
    end
end

for i=1:node3
    for j = 1:node3
        T53(i,j) = T5(4*(i-1) + 1,4*(j-1) + 1);
    end
end

for i=1:node4
    for j = 1:node4
        T54(i,j) = T5(2*(i-1) + 1,2*(j-1) + 1);
    end
end

clear T1 T2 T3 T4 T5

%calculating errors 
error1 = norm(T11 - T51)/norm(T51);
error2 = norm(T22 - T52)/norm(T52);
error3 = norm(T33 - T53)/norm(T53);
error4 = norm(T44 - T54)/norm(T54);
error5 = norm(T55 - T55)/norm(T55);

dx = [ 0.0125, 0.025, 0.05, 0.1];
error = [error4, error3, error2, error1];

%plotting
figure(6)
loglog(dx,error);
hold on
scatter(dx,error);
xlabel("dx")
ylabel("Error")
title("Log Log graph of error vs dx")
legend("Slope is 2.1404")


p = polyfit(log(dx),log(error),1) ;
slope = p(1);