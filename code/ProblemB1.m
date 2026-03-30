clc
clear
g = @circleg;

H = [1/2 1/4 1/8 1/16 1/32];

energyNorm = zeros(length(H),1);
EnE = zeros(length(H),1);

for i = 1:length(H)
    hm = H(i);
        
    [p, e, t] = initmesh(g, 'hmax', hm);
    g1 = u(p(1, e(1,:)), p(2, e(1,:)));
    A = StiffnessAssembler2D(p, t, @(x, y) 1);
    b = LoadAssembler2D(p, t, @f);
    I = eye(length(p));
    
    A(e(1,:),:) = I(e(1,:),:);
    b(e(1,:))= g1;
  
    uh = A\b;

    err = u(p(1,:), p(2,:)) - uh';
    EnE(i) = sqrt(err * A * err');
    
    if i == 1
        figure;
        pdeplot(p, e, t, 'XYData', uh, 'ZData', uh);
    elseif i == length(H)
        figure;
        pdeplot(p, e, t, 'XYData', uh, 'ZData', uh);
    end
end

gamma = zeros((length(H)-1),1);

for i=1:length(energyNorm)-1
    gamma(i) = log(EnE(i+1)/EnE(i))/log(H(i+1)/H(i));
    disp(gamma(i));
end

figure;
loglog(H, EnE, 'r--', 'LineWidth', 1.5);
hold on;
loglog(H, H.^mean(gamma), 'b-', 'LineWidth', 1.5);
xlabel('h_{max}');
grid on;
legend('EnE', ['h_{max}^{\gamma}, \gamma = ' num2str(gamma(end))], 'Location', 'southeast');
set(gca, 'FontSize', 14);
set(legend, 'FontSize', 14);
hold off;

function z = f(x, y)
    z = 8 * pi^2 * sin(2 * pi * x) .* sin(2 * pi * y);
end

function z = u(x, y)
    z = sin(2 * pi * x) .* sin(2 * pi * y);
end
