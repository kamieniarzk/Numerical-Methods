x(1:22) = 2.^(0:21)/300;
%N = 3
I = eye(3);
for k = 1:22
    A3(:, :, k) = makeMatrix(3, x(k));
    A3invLU(:, :, k) = inverseLU(A3(:, :, k));
    A3invLLT(:, :, k) = inverseLLT(A3(:, :, k));
    
    RMS3LU(k) = RMS((A3(:, :, k)*A3invLU(:, :, k) - I));
    RMS3LUnorm(k) = norm(A3(:, :, k)*A3invLU(:, :, k) - I);
    
    RMS3LLT(k) = RMS((A3(:, :, k)*A3invLLT(:, :, k) - I));
    RMS3LLTnorm(k) = norm(A3(:, :, k)*A3invLLT(:, :, k) - I);
    
    ME3LU(k) = maximumError((A3(:, :, k)*A3invLU(:, :, k) - I));
    ME3LUnorm(k) = norm((A3(:, :, k)*A3invLU(:, :, k) - I), inf);
    ME3LLT(k) = maximumError((A3(:, :, k)*A3invLLT(:, :, k) - I));
    ME3LLTnorm(k) = norm((A3(:, :, k)*A3invLLT(:, :, k) - I), inf);
end

%N = 10
I = eye(10);
for k = 1:22
    A10(:, :, k) = makeMatrix(10, x(k));
    A10invLU(:, :, k) = inverseLU(A10(:, :, k));
    A10invLLT(:, :, k) = inverseLLT(A10(:, :, k));
    RMS10LU(k) = RMS((A10(:, :, k)*A10invLU(:, :, k) - I));
    RMS10LUnorm(k) = norm(A10(:, :, k)*A10invLU(:, :, k) - I);
    
    RMS10LLT(k) = RMS((A10(:, :, k)*A10invLLT(:, :, k) - I));
    RMS10LLTnorm(k) = norm(A10(:, :, k)*A10invLLT(:, :, k) - I);
    
    ME10LU(k) = maximumError((A10(:, :, k)*A10invLU(:, :, k) - I));
    ME10LUnorm(k) = norm((A10(:, :, k)*A10invLU(:, :, k) - I), inf);
    ME10LLT(k) = maximumError((A10(:, :, k)*A10invLLT(:, :, k) - I));
    ME10LLTnorm(k) = norm((A10(:, :, k)*A10invLLT(:, :, k) - I), inf);

end

%N = 20
I = eye(20);
for k = 1:22
    A20(:, :, k) = makeMatrix(20, x(k));
    A20invLU(:, :, k) = inverseLU(A20(:, :, k));
    A20invLLT(:, :, k) = inverseLLT(A20(:, :, k));
    RMS20LU(k) = RMS((A20(:, :, k)*A20invLU(:, :, k) - I));
    RMS20LUnorm(k) = norm(A20(:, :, k)*A20invLU(:, :, k) - I);
    
    RMS20LLT(k) = RMS((A20(:, :, k)*A20invLLT(:, :, k) - I));
    RMS20LLTnorm(k) = norm(A20(:, :, k)*A20invLLT(:, :, k) - I);
    
    ME20LU(k) = maximumError((A20(:, :, k)*A20invLU(:, :, k) - I));
    ME20LUnorm(k) = norm((A20(:, :, k)*A20invLU(:, :, k) - I), inf);
    ME20LLT(k) = maximumError((A20(:, :, k)*A20invLLT(:, :, k) - I));
    ME20LLTnorm(k) = norm((A20(:, :, k)*A20invLLT(:, :, k) - I), inf);

end

figure
semilogy(x, RMS3LU, x, RMS10LU, x, RMS20LU)
hold on
title('Root-mean-square for LU')
xlabel('x')
ylabel('RMS')
legend('N=3', 'N=10', 'N=20')

figure
semilogy(x, RMS3LLT, x, RMS10LLT, x, RMS20LLT)
hold on
title('Root-mean-square for LLT')
xlabel('x')
ylabel('RMS')
legend('N=3', 'N=10', 'N=20')

figure
semilogy(x, ME3LU, x, ME10LU, x, ME20LU)
hold on
title('Maximum error for LU')
xlabel('x')
ylabel('Maximum error')
legend('N=3', 'N=10', 'N=20')

figure
semilogy(x, ME3LLT, x, ME10LLT, x, ME20LLT)
hold on
title('Maximum error for LLT')
xlabel('x')
ylabel('Maximum error')
legend('N=3', 'N=10', 'N=20')

% max(A3invLU(:, : , 1)-inv(A3(:, :, 1)))
% max(A3invLU(:, : , 10)-inv(A3(:, :, 10)))
% max(A3invLU(:, : , 20)-inv(A3(:, :, 20)))
% 
% max(A10invLU(:, : , 1)-inv(A10(:, :, 1)))
% max(A10invLU(:, : , 10)-inv(A10(:, :, 10)))
% max(A10invLU(:, : , 20)-inv(A10(:, :, 20)))
% 
% max(A20invLU(:, : , 1)-inv(A20(:, :, 1)))
% max(A20invLU(:, : , 10)-inv(A20(:, :, 10)))
% max(A20invLU(:, : , 20)-inv(A20(:, :, 20)))
% 
% max(A3invLLT(:, : , 1)-inv(A3(:, :, 1)))
% max(A3invLLT(:, : , 10)-inv(A3(:, :, 10)))
% max(A3invLLT(:, : , 20)-inv(A3(:, :, 20)))
% 
% max(A10invLLT(:, : , 1)-inv(A10(:, :, 1)))
% max(A10invLLT(:, : , 10)-inv(A10(:, :, 10)))
% max(A10invLLT(:, : , 20)-inv(A10(:, :, 20)))
% 
% max(A20invLLT(:, : , 1)-inv(A20(:, :, 1)))
% max(A20invLLT(:, : , 10)-inv(A20(:, :, 10)))
% max(A20invLLT(:, : , 20)-inv(A20(:, :, 20)))