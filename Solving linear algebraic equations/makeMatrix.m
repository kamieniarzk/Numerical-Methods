function [Y] = makeMatrix(N, x)
Y(1,1) = x^2;
for i = 2:N
    %filling first column
    Y(i,1) = (-1)^i*3*x/2;
    %filling first row
    Y(1,i) = (-1)^i*3*x/2;
end


for a = 2:N
    for b = 2:N
        %filling rows 
        Y(a, b) = b*9/(4*(-1)^(b-a));
        %inverting rows and columns
        Y(b, a) = Y(a, b);
    end
end


end

