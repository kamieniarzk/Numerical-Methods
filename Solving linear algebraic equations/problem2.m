alfa = linspace(0.99, 1.01, 1000);
    
for i = 1:1000
    y3(1,i) = cond(makeMatrix(3, log10(alfa(1,i))));
    y10(1,i) = cond(makeMatrix(10, log10(alfa(1,i))));
    y20(1,i) = cond(makeMatrix(20, log10(alfa(1,i))));
    
    z3(1,i) = det(makeMatrix(3, log10(alfa(1,i))));
    z10(1,i) = det(makeMatrix(10, log10(alfa(1,i))));
    z20(1,i) = det(makeMatrix(20, log10(alfa(1,i))));
end
figure
semilogy(alfa, y3, alfa, y10, alfa, y20)
xlabel('alfa')
ylabel('cond')
title('Dependence of cond() on alfa')
legend('N=3', 'N=10', 'N=20')
hold on


figure
semilogy(alfa, z3, alfa, z10, alfa, z20)
xlabel('alfa')
ylabel('determinant')
title('Dependence of det() on alfa')
legend('N=3', 'N=10', 'N=20')
hold on