clc 
clear

%number of elements in the domain 
nx = 100; 

%total nodes 
nodes = nx +1; 

%Length of the domain 
Length  = 1; 

%dx 
dx = Length/nx;

%time step 
dt = 1e-4; 

%time limit 
TIME = 25;

%Tollerance to find steady state
tollSteady = 1e-10;

%max iterations
maxIter = 1e10;

%initiallising solution vectors 
T2 = zeros(nodes,1); %Solution vector for case 2 
T2exact = zeros(nodes,1); %exact solution for case 2 

%Coordinates of each node 
pos = 0;
for i = 1:nodes
    crd(i) = pos;
    pos = pos + dx;
end
clear pos;

%exact solution at steady state 
T2exact = (crd.*Length.*Length).^(1/3);

%initial condition 
T2 = sin((pi.*crd)./Length) + crd;

%initialising previous vector
T2Prev = T2;

%initiallising information about the setup 
rho= 1200;
c = 0.013;

%initiallising discritising variable used in the expression 
A = (rho*c*dx*dx)/dt;

%time loop starts for case 2

for i = 1:maxIter
    time = i * dt;
    % space loop starts here
    for j = 2:nodes-1
        T2(j) = (A - 2*T2(j)^2) * T2(j) + (T2(j+1)^2 - T2(j-1)^2)*0.25*(T2(j+1) - T2(j-1)) + T2(j)^2*T2(j-1) + T2(j)^2*T2(j+1);
        T2(j) = T2(j)/A;
    end
    comp = abs((norm(T2) - norm(T2Prev))/(norm(T2) + 1e-12));
%     if(comp < tollSteady)
%         break
%     end
    if(time >= TIME)
        break
    end
    T2Prev = T2;
end

%error with exact solution for case 2 
Error1 = norm(T2 -T2exact) / norm(T2exact);

figure(2)
plot(crd, T2, 'DisplayName' , 'Numerical Solution'  );
hold on 
scatter(crd,T2exact, 'DisplayName', 'Exact Steady Solution');
xlim tight
ylim tight
legend
xlabel("x")
ylabel("Temprature")
title("FTCS - Case 2 - at Time = 25s")




















