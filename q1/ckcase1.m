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

%Tollerance for non linear and gauss sidel loops 
toll = 1e-10;

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

%initialising K vectors 
k = T1;
kPrev = k;

%Time loop starts for case 1 
for i = 1:maxIter
        time = i * dt;
        T1Prev = T1;
        %non linear loop  
        for l = 1:maxIter
                Tk = T1;
                k = Tk;
                kPrev = T1Prev;
                %Gauss sidel loop 
                for n = 1:maxIter
                        Tg = T1;
                        % space loop starts here
                        for j = 2 : nodes - 1
                            T1(j) = A*T1Prev(j) + 0.125*(k(j+1) -  k(j-1))*(T1(j+1) - T1(j-1)) + 0.5*k(j)*(T1(j+1) - 2*T1(j) + T1(j-1)) +...
                                0.125*(kPrev(j+1) - kPrev(j-1))*(T1Prev(j+1) - T1Prev(j-1)) + 0.5*kPrev(j)*(T1Prev(j+1) - 2*T1Prev(j) + T1Prev(j-1));
                            T1(j) = T1(j)/A;
                        end 
                        errorGs = abs(norm(T1-Tg)/norm(T1));
                        if(errorGs < toll)
                            break
                        end
                end
                errorNl = abs(norm(T1-Tk)/norm(T1));
                if(errorNl < toll)
                    break
                end
        end
        error = abs(norm(T1-T1Prev)/norm(T1));
        if(time>=TIME)
            break
        end
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
title("Crank Nicolson - Case 1 - at time = 25s")



