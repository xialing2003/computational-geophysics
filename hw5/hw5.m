clear;
f = 0; T1 = 1; q0 = 1;

sample = 101;
x = linspace(0, 1, sample);
T = zeros(1, sample);

%% get the analytic solution
for i = 1:sample
    T(i) = T1 + (1-x(i))*q0 + 0.5*(1-x(i)^2)*f;
end

figure
plot(x,T,'k');
hold on

%% choose Ne = 2
Ne = 2; 
dx = 1/Ne;
K = zeros(Ne, Ne);
for i = 1:Ne-1
    K(i,i) = K(i,i) + 1/dx;
    K(i,i+1) = -1/dx;
    K(i+1, i) = -1/dx;
    K(i+1,i+1) = K(i+1, i+1) + 1/dx;
end
K(Ne, Ne) = K(Ne, Ne) + 1/dx;

F = ones(Ne, 1)*dx*f;
F(1) = 0.5*dx*f + q0; 
F(Ne) =  dx*f + T1/dx; 

dB = inv(K)*F;
dB(Ne+1) = T1;

for i = 1:sample-1
    cx = x(i); %current x
    e = floor(cx/dx) + 1;
    % disp(e);
    T(i) = dB(e)*(e*dx - cx)/dx + dB(e+1)*(cx - (e-1)*dx)/dx;
end
T(101) = T1;
plot(x,T,'--r');
hold on

%% choose Ne = 10

Ne = 10; 
dx = 1/Ne;
K = zeros(Ne, Ne);
for i = 1:Ne-1
    K(i,i) = K(i,i) + 1/dx;
    K(i,i+1) = -1/dx;
    K(i+1, i) = -1/dx;
    K(i+1,i+1) = K(i+1, i+1) + 1/dx;
end
K(Ne, Ne) = K(Ne, Ne) + 1/dx;

F = ones(Ne, 1)*dx*f;
F(1) = 0.5*dx*f + q0; 
F(Ne) =  dx*f + T1/dx; 

dB = inv(K)*F;
dB(Ne+1) = T1;

for i = 1:sample-1
    cx = x(i); %current x
    e = floor(cx/dx) + 1;
    % disp(e);
    T(i) = dB(e)*(e*dx - cx)/dx + dB(e+1)*(cx - (e-1)*dx)/dx;
end
T(101) = T1;
plot(x,T,'--g');

xlabel('location x')
ylabel('temperature T')
title('heat equation colution when f = 0, T_1 = 1, q_0 = 1')
legend('analytic', 'Ne=2', 'Ne=10');
saveas(gcf, 'f=0.png');

