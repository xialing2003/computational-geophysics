clear;

%% exercise 1 
L = pi/2;f = 0; TL = 1; kapa = 1; rhocp = 1; q0 = 0;

% construct matrices and the force
Ne = 10; 
dx = L/Ne;
K = zeros(Ne, Ne);
for i = 1:Ne-1
    K(i,i) = K(i,i) + kapa/dx;
    K(i,i+1) = -kapa/dx;
    K(i+1, i) = -kapa/dx;
    K(i+1,i+1) = K(i+1, i+1) + kapa/dx;
end
K(Ne, Ne) = K(Ne, Ne) + kapa/dx;

M = zeros(Ne, Ne);
for i = 1:Ne-1
    M(i,i) = M(i,i) + dx*rhocp/3;
    M(i,i+1) = dx*rhocp/6;
    M(i+1, i) = dx*rhocp/6;
    M(i+1,i+1) = M(i+1, i+1) + dx*rhocp/3;
end
M(Ne, Ne) = M(Ne, Ne) + dx*rhocp/3;

F = ones(Ne, 1)*dx*f;
F(1) = 0.5*dx*f + q0; 
F(Ne) =  dx*f + kapa/dx*TL; 

% get the initial coefficients
d = zeros(Ne+1,1);
d(1) = 2;
for i = 2:Ne
    point_x = (i-1)*dx;
    d(i) = 1 + cos(point_x);
end
d(Ne + 1) = TL;
v = M\(F-K*d(1:Ne));

% temperature changes with time
sample = 101;
x = linspace(0, L, sample);
T_simu = zeros(1, sample);
T_ana = zeros(1, sample);

for i = 1:sample
    T_ana(i) = 1 + cos(x(i));
end
for i = 1:sample-1
    e = floor(x(i)/dx) + 1;
    T_simu(i) = d(e)*(e*dx - x(i))/dx + d(e+1)*(x(i) - (e-1)*dx)/dx;
end
T_simu(sample) = TL;

figure('Position',[100 100 500 600])
subplot(6,1,1)
plot(x, T_ana, 'k')
hold on
plot(x, T_simu,'--r')
title('time = 0s')
ylabel('Temp')

dt = 0.1; Nt = 50;
for k = 1:Nt
    d(1:Ne) = d(1:Ne) + 0.5 * dt * v;
    v = (M+0.5*dt*K)\(F - K*d(1:Ne));
    d(1:Ne) = d(1:Ne) + 0.5 * dt * v;
    
    for i = 1:sample
        T_ana(i) = 1 + exp(-dt*k)*cos(x(i));
    end
    for i = 1:sample-1
        e = floor(x(i)/dx) + 1;
        T_simu(i) = d(e)*(e*dx - x(i))/dx + d(e+1)*(x(i) - (e-1)*dx)/dx;
    end
    T_simu(sample) = TL;
    if mod(k,10) == 0
        subplot(6, 1, k/10+1);
        plot(x, T_ana, 'k')
        hold on
        plot(x, T_simu,'--r')
        title(['time = ',num2str(k*dt),'s'])
        ylabel('Temp')
    end
end
xlabel('location x')
saveas(gcf, 'exercise1.png');

%% exercise 2 

L = 20;f = 0; TL = 0; kapa = 1; rhocp = 1; q0 = 0;
% construct matrices and the force
Ne = 100; 
dx = L/Ne;
K = zeros(Ne, Ne);
for i = 1:Ne-1
    K(i,i) = K(i,i) + kapa/dx;
    K(i,i+1) = -kapa/dx;
    K(i+1, i) = -kapa/dx;
    K(i+1,i+1) = K(i+1, i+1) + kapa/dx;
end
K(Ne, Ne) = K(Ne, Ne) + kapa/dx;

M = zeros(Ne, Ne);
for i = 1:Ne-1
    M(i,i) = M(i,i) + dx*rhocp/3;
    M(i,i+1) = dx*rhocp/6;
    M(i+1, i) = dx*rhocp/6;
    M(i+1,i+1) = M(i+1, i+1) + dx*rhocp/3;
end
M(Ne, Ne) = M(Ne, Ne) + dx*rhocp/3;

F = ones(Ne, 1)*dx*f;
F(1) = 0.5*dx*f + q0; 
F(Ne) =  dx*f + kapa/dx*TL; 

% get the initial coefficients
d = ones(Ne+1,1);
d(Ne + 1) = TL;
v = M\(F-K*d(1:Ne));

% temperature changes with time
sample = 101;
x = linspace(0, L, sample);
T_simu = zeros(1, sample);
T_ana = zeros(1, sample);

for i = 1:sample
    T_ana(i) = 1;
end
for i = 1:sample-1
    e = floor(x(i)/dx) + 1;
    T_simu(i) = d(e)*(e*dx - x(i))/dx + d(e+1)*(x(i) - (e-1)*dx)/dx;
end
T_simu(sample) = TL;

figure('Position',[100 100 500 900])
subplot(10,1,1)
plot(x, T_ana, 'k')
hold on
plot(L-x, T_simu,'--r')
title('time = 0s')
ylabel('Temp')
%
dt = 0.1; Nt = 90;
for k = 1:Nt
    d(1:Ne) = d(1:Ne) + 0.5 * dt * v;
    v = (M+0.5*dt*K)\(F - K*d(1:Ne));
    d(1:Ne) = d(1:Ne) + 0.5 * dt * v;
    
    for i = 1:sample
        value = x(i)/(2*sqrt(kapa/rhocp*dt*k));
        T_ana(i) = 1 - erfc(value);
    end
    for i = 1:sample-1
        e = floor(x(i)/dx) + 1;
        T_simu(i) = d(e)*(e*dx - x(i))/dx + d(e+1)*(x(i) - (e-1)*dx)/dx;
    end
    T_simu(sample) = TL;
    if mod(k,10) == 0
        subplot(10, 1, k/10+1);
        plot(x, T_ana, 'k')
        hold on
        plot(L-x, T_simu,'--r')
        title(['time = ',num2str(k*dt),'s'])
        ylabel('Temp')
    end
end
xlabel('location x')
saveas(gcf, 'exercise2.png');