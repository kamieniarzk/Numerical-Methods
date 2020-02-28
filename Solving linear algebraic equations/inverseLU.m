function [X] = inverseLU(A)

[L, U, P] = lu(A);
n = size(A, 1);
I = eye(n);
%LY = I, therefore Y = I/L
for c = 1:n %column number
    for r = 1:n %row number
        sum = 0;
        for k = 1:(r-1) %for summation 
            sum = sum + L(r, k)*Y(k, c);
        end
        Y(r, c) = (I(r, c) - sum)/L(r, r);
    end
    
    for r = 1:n
        sum = 0;
        j = r - 1;
        for k = (n - j + 1):n %to go in reverse
            sum = sum + U(n - j, k)*X(k, c);
        end
        X(n - j, c) = (Y(n - j, c) - sum)/U(n - j, n - j);
    end
end

X = X*P;

