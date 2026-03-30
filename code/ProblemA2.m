a = -1;
b = 1;
gl = 0;
gr = 0;
alpha = 0.01;

beta = 0.9;
N = 11;
TOL = 1e-3;
Nmax = 1e4;
h = (b - a) / N;
x = a:h:b;

while true
    N = length(x);
    zeta = my_discrete_Lapl(x, @fStar, gl, gr);
    eta2 = zeros(N, 1);
    for j = 1 : N-1
        h = x(j+1)- x(j);
        a = f(x(j)) + alpha*zeta(j);
        b = f(x(j+1)) + alpha*zeta(j+1);
        t = (a.^2+b.^2)*h/2;
        eta2(j) = h.^2 * t;
    end
    
    if sum(eta2) <= TOL || N > Nmax
        break;
    end
    
    for k = 1 : (length(eta2)-1)
        if eta2(k) > beta*max(eta2)
            x = [x (x(k+1)+x(k))/2];
        end
    end
    x = sort(x);
end

A = my_stiffness_matrix_assemble(x);
B = my_load_vector_assemble(x, @fStar, gl, gr);
xi = A \ B;
Re = f(x) + alpha * zeta';
eta = eta2.^0.5;

figure;
subplot(2, 2, 1);
plot(x, xi);
title('Solution uh');
xlabel('x');
ylabel('uh');

% Plot the residual R(uh)
subplot(2, 2, 2);
plot(x, Re);
title('Residual R(uh)');
xlabel('x');
ylabel('R(uh)');

% Plot the error indicator eta(uh)
subplot(2, 2, 3);
plot(x, eta);
title('Error Indicator eta(uh)');
xlabel('x');
ylabel('eta(uh)');

% Plot the mesh-size distribution
subplot(2, 2, 4);
plot(x(2:end), 1./diff(x));
title('Mesh-size Distribution');
xlabel('x');
ylabel('1/dx');

disp(size(x));


function y = f(x)
    condition = ((x <= 0.8) & (x >= 0.2)) | ((x <= -0.2) & (x >= -0.8));
    y = zeros(size(x));
    y(condition) = 10;
end

function y = fStar(x)
    condition = ((x <= 0.8) & (x >= 0.2)) | ((x <= -0.2) & (x >= -0.8));
    y = zeros(size(x));
    y(condition) = 10 / 0.01;
end

function zeta = my_discrete_Lapl(x, ~, gl, gr)
    A = my_stiffness_matrix_assemble(x);
    B = my_load_vector_assemble(x, @fStar, gl, gr);
    xi = A \ B;
    M = my_mass_matrix_assemble(x);
    zeta = -inv(M) * A * xi; 
end

function A = my_stiffness_matrix_assemble(x)

    N = length(x) - 1;
    A = zeros(N+1, N+1); 
    
    for i = 1 : N
        h = x(i+1) - x(i);
        n = [i i+1];
        A(n,n) = A(n,n) + [1 -1; -1 1]/h;
    end
    
    % adjust for BC
    A(1,1) = 1;
    A(1,2) = 0;
    A(N+1,N) = 0;
    A(N+1,N+1) = 1;

end

function b = my_load_vector_assemble(x, f, gl, gr)
    %
    % Return assembled load vector b
    % Input vector x of node coords
    % 
    N = length(x) - 1;
    b = zeros(N+1, 1);
    for i = 1 : N
        h = x(i+1) - x(i);
        n = [i i+1];
        b(n) = b(n) + [f(x(i)); f(x(i+1))]*h/2;
    end
    b(1) = gl;
    b(N+1) = gr;

end

function M = my_mass_matrix_assemble(x)
    N = length(x)-1;
    M = zeros(N+1, N+1);
    
    for i = 1 : N
        h = x(i+1) - x(i);
        n = [i i+1];
        M(n,n) = M(n,n) + [1/3 1/6; 1/6 1/3]*h;
    end
end
