clc 
clear

%dx value given 
dx = 0.1;

%dy value given 
dy = 0.1;

%Length given - 
Length = 1;

%No of elements - 
elemX = Length/dx; %x axis
elemY = Length/dy; %y axis

%Number of nodes
nodesX = elemX + 1;
nodesY = elemY + 1;

%information given 
alpha = 0.0001;

%dt 
dt = 0.01;

%time to run the simulation 
TIME = 10;

%Max iterations limit 
maxIter = 1e10;

%Tollerance 
toll = 1e-6;

%initiallising discritising variable used in the expression 
A = (dx*dx)/(dt*alpha);

%meshgrid formation 
crdx = (1:nodesX);
crdy = (1:nodesY);
crdx = (crdx-1)/elemX;
crdy = (crdy-1)/elemY;

%initial conditions
for(i = 1:nodesX)
    for(j=1:nodesY)
        if(j>(nodesY+1)/2 && i>(nodesY+1)/2)
            continue
        end
        T(i,j) = 400*sin(pi*crdx(i))*sin(pi*crdy(j)) + 400*crdx(i)*crdy(j);
    end
end

%Initiating neumann boundary conditions
for j=1:(nodesY+1)/2
    T(nodesX+1,j) = T(nodesX-1,j);
end
for j=1:(nodesX+1)/2
    T(j,nodesY+1) = T(j,nodesY-1);
end

%Initiallising temp initial vector 
Tinitial = T;

%temporal loop starts here 
for t = 1:maxIter
    time = t*dt;
    TPrev = T;
    %Gauss Sidel loop starts here 
    for l = 1:maxIter
        Tg = T;
        %space loop starts here 
        for m = 2:nodesX
            for n = 2:nodesY
                if(n>=(nodesY+1)/2 && m>=(nodesY+1)/2)
                    continue
                end
                T(m,n) = TPrev(m,n)*A + 0.5*(TPrev(m+1,n) - 2*TPrev(m,n) + TPrev(m-1,n)) +...
                    0.5*(TPrev(m,n+1) - 2*TPrev(m,n) + TPrev(m,n-1)) +...
                    0.5*(T(m+1,n) - 2*T(m,n) + T(m-1,n)) + 0.5*(T(m,n+1) - 2*T(m,n) + T(m,n-1));
                T(m,n) = T(m,n)/A;
            end 
        end
        %Initiating neumann boundary conditions
        for j=1:(nodesY+1)/2
            T(nodesX+1,j) = T(nodesX-1,j);
        end
        for j=1:(nodesX+1)/2
            T(j,nodesY+1) = T(j,nodesY-1);
        end
        errorGs = abs(norm(T - Tg)/(norm(T)+1e-12));
        if(errorGs < toll)
            break
        end
    end
    errorSt = abs(norm(T - TPrev)/(norm(T)+1e-12));
    if(time >= TIME)
        break 
    end
end

%Plotting 
for i = 1:(nodesY+1)/2
    pltT(i) = T((nodesX+1)/2 , i);
    pltY(i) = crdy(i);
end


figure(1)
plot(pltY,pltT ,'DisplayName' , 'Numerical Solution');
hold on
scatter(pltY,pltT, 'DisplayName' , 'Data points');
hold off
legend
xlabel("x")
ylabel("Temprature")
title("dx = 0.1 - at time = 10s")

T1 = T;
save x01.mat T1
















