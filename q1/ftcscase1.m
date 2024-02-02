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
T1 = zeros(nodes,1); %solution vector for case 1 
T1exact = zeros(nodes,1); %exact solution for case 1 

%Coordinates of each node 
pos = 0;
for i = 1:nodes
    crd(i) = pos;
    pos = pos + dx;
end
clear pos;

%exact solution at steady state 
T1exact = sqrt(crd.*Length);

%initial condition 
T1 = sin((pi.*crd)./Length) + crd;

%initialising previous vector
T1Prev = T1;

%initiallising information about the setup 
rho= 1200;
c = 0.013;

%initiallising discritising variable used in the expression 
A = (rho*c*dx*dx)/dt;

%time loop starts for case 1

for i = 1:maxIter
    time = i * dt;
    % space loop starts here
    for j = 2:nodes-1
        T1(j) = (A - 2*T1(j)) * T1(j) + (T1(j+1) - T1(j-1))*0.25*(T1(j+1) - T1(j-1)) + T1(j)*T1(j-1) + T1(j)*T1(j+1);
        T1(j) = T1(j)/A;
    end
    comp = abs((norm(T1) - norm(T1Prev))/(norm(T1) + 1e-12));
%     if(comp < tollSteady)
%         break
%     end
    if time >= TIME
        break 
    end
    T1Prev = T1;
end

%error with exact solution for case 1 
Error1 = norm(T1 -T1exact) / norm(T1exact);

figure(1)
plot(crd, T1, 'DisplayName' , 'Numerical Solution'  );
hold on 
scatter(crd,T1exact, 'DisplayName', 'Exact Steady Solution');
xlim tight
ylim tight
legend
xlabel("x")
ylabel("Temprature")
title("FTCS - Case 1 - at time = 25s")





















