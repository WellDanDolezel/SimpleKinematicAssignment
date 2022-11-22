clear, clc, close all

a = 0.1;
b = 0.2;
fi_0 = pi/6;
t = linspace(0, 1, 101);
omega = - 1;
init_estim = [0; 0];
tol = 1e-6;

%% Solving positions
for ii = 1:length(t)
    fi = fi_0 + omega * t(ii);

    F_pos = @(x) [a*cos(fi) + b*cos(x(1)) - x(2);
                  a*sin(fi) - b*sin(x(1))];

    J_pos = @(x) [- b*sin(x(1)), -1;
                  - b*cos(x(1)), 0];

    [x(:,ii),n] = NR_function(F_pos, J_pos, init_estim, tol);
    init_estim = x(:,ii);
end


theta = x(1,:);
d = x(2,:);

figure(1), set(gcf,'color','w','Position',[200 200 600 175]);
colororder({'b','r'})
yyaxis left
plot(t, theta, 'LineWidth', 2)
grid on, grid minor
xlabel ('t [s]','Interpreter','latex')
ylabel ('${\theta}$ [rad]','Interpreter','latex')

yyaxis right
plot(t, d, 'LineWidth', 2)
ylabel ('$d$ [m]','Interpreter','latex')



%% Solving velocities

init_estim = [0; 0];
for ii = 1:length(t)
    fi = fi_0 + omega * t(ii);

    F_vel = @(x) [- a*sin(fi)*omega - b*sin(theta(ii))*x(1) - x(2);
                    a*cos(fi)*omega - b*cos(theta(ii))*x(1)];

    J_vel = @(x) [- b*sin(theta(ii)), -1;
                  - b*cos(theta(ii)), 0];

    [x(:,ii),n] = NR_function(F_vel, J_vel, init_estim, tol);
    init_estim = x(:,ii);
end

theta_dot = x(1,:);
d_dot = x(2,:);


figure(2), set(gcf,'color','w','Position',[200 200 600 175]);
colororder({'b','r'})
yyaxis left
plot(t, theta_dot, 'LineWidth', 2)
grid on, grid minor
xlabel ('t [s]','Interpreter','latex')
ylabel ('$\dot{\theta}$ [rad/s]','Interpreter','latex')

yyaxis right
plot(t, d_dot, 'r', 'LineWidth', 2)
ylabel ('$\dot{d}$ [m/s]','Interpreter','latex')

