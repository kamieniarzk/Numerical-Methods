%---ASSIGNMENT C ENUME---
%Kacper Kamieniarz
%293065

tspan = 0:0.01:10; %integration step h = 0.01

%---ODE113----------------------------------------------
options = odeset('RelTol', 1e-16, 'AbsTol', 1e-16);
[t113,y113] = ode113(@f, [0 10],[0 2], options);
figure
plot(t113, y113(:, 1), t113, y113(:, 2))
% title('Solution of differential equation using ODE113');
xlabel('t');
ylabel('y');
% legend('y_1','y_2')

%---MyODE----------------------------------------------

% figure
[t, y] = myODE(tspan);
plot(t, y)
hold on
title('Solution of differential equation using MyODE');
xlabel('t');
ylabel('y');
% legend('y_1','y_2')

%---TASK 2----------------------------------------------
h = logspace(-4, -1, 20);

for i = 1:20
    RG6RMS(i) = RMS(h(i), @myODE);
    RG6ME(i) = ME(h(i), @myODE);
    eulerRMS(i) = RMS(h(i), @euler);
    eulerME(i) = ME(h(i), @euler);
end

figure
loglog(h, RG6RMS,'-og')
hold on 
loglog(h, eulerRMS,  '-ob')
title('Dependency of Root-Mean-Square error on step h');
xlabel('h');
ylabel('Error');
legend('Runge-Kutta 6','Euler explicit')

figure
loglog(h, RG6ME, '-og')
hold on 
loglog(h, eulerME,  '-ob')
title('Dependency of Maximum Error on step h');
xlabel('h');
ylabel('Error');
legend('Runge-Kutta 6','Euler explicit')


% 
% % 
%---EULER----------------------------------------------
[te, ye] = euler(tspan);
figure
plot(te, ye)
title('Solution of differential equation using Explicit Forward Eulers method');
xlabel('t');
ylabel('y');
legend('y_1','y_2')

%---FUNCTIONS-------------------------------------------

%---f(t,y)---
function dydt = f(t, y)
    dydt = [y(2); -2/3*y(2)-10/9*y(1)];
    %basing on y' = A*y
end

%---forward Euler explicit---
function [x, u] = euler(tspan)
    x = tspan;
    u = [0; 2]; %Initial conditions
    h = tspan(2) - tspan(1);
    for i = 1:length(tspan) - 1
        u(:, i+1) = u(:, i) + h*f(x(i), u(:, i));
    end
end
    
%---Runge-Kutta order 6 Gauss-Legendre---
function [x, u] = myODE(tspan)
        x = tspan;
        A = [0 1; -10/9 -2/3]; %from y' = A*y
        u = [0; 2]; %initial conditions
        h = tspan(2) - tspan(1); %integration step
        I = eye(2);
        M = [I - 5/36*h*A, -(2/9 - sqrt(15)/15)*h*A, -(5/36 - sqrt(15)/30)*h*A;
        -(5/36 + sqrt(15)/24)*h*A, I - 2/9*h*A, -(5/36 - sqrt(15)/24)*h*A;
        -(5/36 + sqrt(15)/30)*h*A, -(2/9 + sqrt(15)/15)*h*A, I - 5/36*h*A];
        for i = 1:length(tspan) - 1
            R = [A*u(:, i); A*u(:, i); A*u(:, i)];
            k = M\R;
            g = [k(1) k(3) k(5); k(2) k(4) k(6)]; %k1, k2, k3 coefficients
            u(:, i+1) = u(:, i) + h*(5/18*g(:, 1) + 4/9*g(:, 2) +5/18*g(:, 3));
        end
end
            
%---RMS---
function error = RMS(h, func)
    %Setting the tolerance parameters
    options = odeset('RelTol', 1e-16, 'AbsTol', 1e-16);
    span = 0:h:10;
    [t113, y113] = ode113(@f, span,[0 2], options);
    y113 = y113'; %transpose to column vector
    [x, u] = func(span);
    error = norm((u(1, :) - y113(1, :)), 2)/norm(y113(1, :), 2);
end

%---ME---
function error = ME(h, func)
    %Setting the tolerance parameters
    options = odeset('RelTol', 1e-16, 'AbsTol', 1e-16);
    span = 0:h:10;
    [t113, y113] = ode113(@f, span,[0 2], options);
    y113 = y113'; %transpose to column vector
    [x, u] = func(span);
    error = norm((u(1, :) - y113(1, :)), 'Inf')/norm(y113(1, :), 'Inf');
end