function [B] = RMS(A)
X = A*A.';
B = max(sqrt(eig(X)));
end
