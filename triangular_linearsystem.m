function U = triangular_system(A,f);
% Function U = trinagular_system is used to solve a linear system with a triangular matrix A and 
% right hand site vector f.

M = size(A,1);
epsilon = zeros(1,M-1);
g = epsilon;
epsilon(1) = - A(1,2)/A(1,1);
g(1) = f(1)/A(1,1);
% calculate intermediate values: epsilon and g.
for k=2:M-1;
    var = A(k,k-1)*epsilon(k-1) + A(k,k);
    epsilon(k) = -A(k,k+1)/var;
    g(k)  = (f(k) - A(k,k-1)*g(k-1))/var;
end
% Calculate the solution.
solution = zeros(M,1);
solution(M) = (f(M) - A(M,M-1)*g(M-1))/(A(M,M-1)*epsilon(M-1) + A(M,M));
for k=M-1:-1:1
    solution(k) = epsilon(k)* solution(k+1) + g(k);
end
U = solution;
