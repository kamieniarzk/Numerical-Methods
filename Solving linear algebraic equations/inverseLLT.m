function [X] = inverseLLT(A)

n = size(A, 1);
I = eye(n);
L = chol(A);
LT = transpose(L);
for c = 1:n %column number
    for r = 1:n %row number
        sum = 0;
        for k = 1:(r-1) %for summation 
            sum = sum + LT(r, k)*Y(k, c);
        end
        Y(r, c) = (I(r, c) - sum)/LT(r, r);
    end
    
    for r = 1:n
        sum = 0;
        j = r - 1;
        for k = (n - j + 1):n %to go in reverse
            sum = sum + L(n - j, k)*X(k, c);
        end
        X(n - j, c) = (Y(n - j, c) - sum)/L(n - j, n - j);
    end
end





end

