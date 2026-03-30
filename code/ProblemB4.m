clc
clear

alpha = 0.01;
beta = 1;
gamma = 0.2;

T = 30; % final time

dt = 0.1; % time step
tn = T / dt; % number of time levels
MLoss = zeros(tn+1, 2);

g = @circleg;

H = [1/5, 1/20];

for index = 1 : length(H)
    hmax = H(index); % mesh size
    [p, e, t] = initmesh(g, 'hmax', hmax);
    nt = size(t, 2); % number of elements
    
    A = StiffnessAssembler2D(p, t, @(x, y) alpha);
    M = MassAssembler2D(p, t);
    I = eye(length(p));
    
    A(e(1,:),:) = I(e(1,:),:);

    old_xi = u0(p(1,:), p(2,:))';
    for i = 1: tn
        D = ((1/dt).*M + (1/2).*A);
        D(e(1,:),:) = I(e(1,:),:);
        b = (((1/dt).*M - (1/2).*A)*old_xi + beta.*M*old_xi - beta*gamma.*M*(old_xi.*old_xi));
        b(e(1,:),:) = 0;
        xi = D \ b;

        % compute the mass loss
        for K = 1:nt
            triangle = t(1:3, K);
            x = p(1, triangle);
            y = p(2, triangle);
            [area, ~, ~] = HatGradients(x, y);
            MLoss(i+1, index) = MLoss(i+1,index) + area/3 * (sum(u0(x,y))- sum(xi(triangle)));
        end
        old_xi = xi;
    end
end

figure;
plot(linspace(0, T, tn + 1), MLoss(:,1), 'DisplayName', 'h_{max} = 1/5', 'LineWidth', 2, 'Color', 'b');
hold on;
plot(linspace(0, T, tn + 1), MLoss(:,2), 'DisplayName', 'h_{max} = 1/20', 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--');
hold off;
legend('show');
xlabel('Time (day)', 'Interpreter', 'latex');
ylabel('Concentration of Hormone ($mmol/mm^3$)', 'Interpreter', 'latex');

figure;
pdeplot(p, e, t, 'XYData', u0(p(1,:), p(2,:)), 'ZData', u0(p(1,:), p(2,:)));

figure;
pdeplot(p, e, t, 'XYData', xi, 'ZData', xi);
