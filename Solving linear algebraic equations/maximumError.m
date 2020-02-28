function [B] = maximumError(A)
B = max((sum(abs(A), 2)));
end

