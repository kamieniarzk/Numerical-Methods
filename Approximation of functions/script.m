hold off;
close all;
clear;
clc;
% TASK 1 
% -------------------------------------------------------
x = linspace(-1, 1, 500)';
f = makeF(x);

% N = 10
x10 = linspace(-1, 1, 10)';
y10 = makeF(x10);

% N = 20
x20 = linspace(-1, 1, 20)';
y20 = makeF(x20);

% N = 30
x30 = linspace(-1, 1, 30)';
y30 = makeF(x30);

% TASK 2 
% -------------------------------------------------------
figure 
plot(x, makeApp(5, 10))
xlabel('x')
ylabel('f(x)')
hold on
plot(x, f)
hold on
plot(x10, y10, 'or')
title('Approx. with K = 5, N = 10')
legend('approximation', 'function', 'yn used for approx')
figure
plot(x, makeApp(9, 10))
xlabel('x')
ylabel('f(x)')
hold on
plot(x, f)
hold on
plot(x10, y10, 'or')
title('Approx. with K = 9, N = 10')
legend('approximation', 'function', 'yn used for approx')
figure 
plot(x, makeApp(5, 20))
xlabel('x')
ylabel('f(x)')
hold on
plot(x, f)
hold on
plot(x20, y20, 'or')
title('Approx. with K = 5, N = 20')
legend('approximation', 'function', 'yn used for approx')
figure
plot(x, makeApp(9, 20))
xlabel('x')
ylabel('f(x)')
hold on
plot(x, f)
hold on
plot(x20, y20, 'or')
title('Approx. with K = 9, N = 20')
legend('approximation', 'function', 'yn used for approx')
figure
plot(x, makeApp(5, 30))
xlabel('x')
ylabel('f(x)')
hold on
plot(x, f)
hold on
plot(x30, y30, 'or')
title('Approx. with K = 5, N = 30')
legend('approximation', 'function', 'yn used for approx')
figure
plot(x, makeApp(9, 30))
xlabel('x')
ylabel('f(x)')
hold on
plot(x, f)
hold on
plot(x30, y30, 'or')
title('Approx. with K = 9, N = 30')
legend('approximation', 'function', 'yn used for approx')

%TASK 3 
%-------------------------------------------------------

for n = 5:50
    for k = 3:n
        error2(k,n) = delta2(k, n);
        errorInf(k,n) = deltaInf(k,n);
    end
end

figure
surf(log10(error2))
title('RMS error for N = 5,...,50 K < N')
xlabel('N')
ylabel('K')
zlabel('Error')
hold on

figure
surf(log10(errorInf))
title('Maximum error error for N = 5,...,50 K < N')
xlabel('N')
ylabel('K')
zlabel('Error')
        
%TASK 4 
% -------------------------------------------------------

power = linspace(-5, -1, 20);
index = 1;
for i = power
    for n = 5:50
        yc = errorData(n, i);
        for k = 3:n-1
            errorc2(k, n) = corrupted2(k, yc); %RMS
            errorcInf(k, n) = corruptedInf(k, yc); %Maximum error
            
            min2(index) = min(errorc2(errorc2 > 0), [], 'all'); %RMS
            minInf(index) = min(errorcInf(errorcInf > 0), [], 'all'); %Maximum error
        end
    end
    index = index +1;
end


deviation = 10.^(power);
%RMS FITTING
g = polyfit(deviation, min2, 5);
gval = polyval(g, deviation);
figure
semilogy(log10(deviation), gval)
hold on
semilogy(log10(deviation), min2, 'or')
hold off
title('Min value of RMS for standard deviation')
xlabel('Sigma power');
ylabel('Minimum of RMS');
legend('polyfit', 'min(RMS)');

%Maximum error fitting
gg = polyfit(deviation, minInf, 2);
ggval = polyval(gg, deviation);
figure
semilogy(log10(deviation), ggval)
hold on
semilogy(log10(deviation), minInf, 'or')
title('Min value of ME for standard deviation')
xlabel('Sigma power');
ylabel('Minimum of Maximum Error');
legend('polyfit', 'min(ME)');


% FUNCTIONS
% --------------------------------------------------------

function [B] = makeF(x) %Approximated function
    B = sqrt(1 - x.^2).*exp(x - 1/3);
end

function [FI] = makeFI(K, N) %Makes FI matrix 
    x = linspace(-1, 1, N)';
    FI(1:N, 1) = 1;
    FI(:, 2) = x;
    if K > 2
        for j = 1:N
                for i = 3:K
                    FI(j, i) = 2*FI(j, i - 1).*x(j) - FI(j, i - 2);
                end
        end
    end
end

function [B] = makeApp(K, N) %Makes approx. for given N and K
    x = linspace(-1, 1, N)';
    y = makeF(x);
    FI = makeFI(K, N);
    p = (FI' * FI)\(FI'*y);
    FI1 = makeFI(K, 500);
    B = FI1 * p;
end

function [B] = delta2(K, N)
    x = linspace(-1, 1, 500)';
    B = norm(makeApp(K, N) - makeF(x))/norm(makeF(x));
end

function [B] = deltaInf(K, N)
    x = linspace(-1, 1, 500)';
    B = norm(makeApp(K, N) - makeF(x), 'inf')/norm(makeF(x), 'inf');
end

function [B] = errorData(N, power)
    %delta y - pseudorandom numbers from normal distribution with variance
    %M^2
    yp = (randn(N, 1)*10^power)';
    x = linspace(-1 ,1, N)';
    B = makeF(x) + yp;
end

function [B] = errorApp(K, yn) 
%Makes approximation but with corrupted data
    N = size(yn, 1);
    x = linspace(-1, 1, N)';
    FI = makeFI(K, N);
    p = (FI' * FI)\(FI'*yn);
    FI1 = makeFI(K, 500);
    B = FI1 * p;
end

function [B] = corrupted2(K, yn) 
%makes delta 2 with corrupted data
    x = linspace(-1, 1, 500)';
    B = norm(errorApp(K, yn) - makeF(x))/norm(makeF(x));
end

function [B] = corruptedInf(K, yn) 
%makes delta inf with corrupted data
    x = linspace(-1, 1, 500)';
    B = norm(errorApp(K, yn) - makeF(x), 'inf')/norm(makeF(x), 'inf');
end



