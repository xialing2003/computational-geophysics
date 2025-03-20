N = 1001;
dx = 0.1;
lc = linspace(0,100,N);
u = zeros(1,N);
for x = 1:N
    u(x) = exp(-0.1*((x-1)*dx-50)^2);
end
kx = fftshift((-N/2:N/2-1) * (2*pi / (N*dx))); 

%% compare the initial stress status calculations 
kapa = 1;
T_fd = zeros(1,N);
T_fd(1) = kapa*(u(2) - u(1))/dx;
T_fd(N) = kapa*(u(N) - u(N-1))/dx;
for x = 2:N-1
    T_fd(x) = kapa*(u(x+1)-u(x-1))/(2*dx);
end

Fu = fft(u);
kFu = 1i * kx .* Fu; 
FT = kapa * ifft(kFu,'symmetric');

figure
plot(lc, T_fd, 'b', 'DisplayName', 'Finite Difference'); 
hold on
plot(lc, real(FT), 'r--', 'DisplayName', 'FFT Derivative'); 
xlabel('location(m)')
title('Two methods to calculate the initial stress status')
legend;
saveas(gcf, 'comp_T.png');

%% homogeneous case
kapa = 1; rho = 1;
v0 = zeros(1,N);
T0 = FT; 
dt = 0.01;

kFT = 1i * kx .* fft(T0);                    % ik F(T)
v1 = v0 + dt * rho * ifft(kFT,'symmetric');   % dv/dt = rho * dT/dx
kFv = 1i * kx .* fft(v0);
T1 = T0 + dt * kapa * ifft(kFv, 'symmetric');

% figure
figure('Units','centimeter','Position',[2 2 25 40]);
subplot(11,2,1)
plot(lc,v0)
ylim([-0.3, 0.3])
title('velocity at time = 0 s');
subplot(11,2,2)
plot(lc,T0)
ylim([-0.3, 0.3])
title('stress at time = 0 s');

for i = 2:10000
    kFT = 1i * kx .* fft(T1);                    % ik F(T)
    v2 = v0 + 2 * dt * rho * ifft(kFT,'symmetric');   % dv/dt = rho * dT/dx
    kFv = 1i * kx .* fft(v1);
    T2 = T0 + 2 * dt * kapa * ifft(kFv, 'symmetric');
    if mod(i,1000) == 0
        subplot(11,2,2*(i/1000)+1)
        plot(lc, real(v2))
        ylim([-0.3, 0.3])
        title(sprintf('velocity at time = %d s', i*dt));
        subplot(11,2,2*(i/1000)+2)
        plot(lc, real(T2))
        ylim([-0.3, 0.3])
        title(sprintf('stress at time = %d s', i*dt));
    end
    v0 = v1; v1 = v2;
    T0 = T1; T1 = T2;
end
xlabel('location(m)')
subplot(11,2,21)
xlabel('location(m)')
saveas(gcf, 'homo.png');

%% inhomogeneous case
kapa1 = 1; kapa2 = 4; rho = 1;
v0 = zeros(1,N);

Fu = fft(u);
kFu = 1i * kx .* Fu; 
dudx = ifft(kFu,'symmetric');
T0(1:60/dx+1) = kapa1 * dudx(1:60/dx+1);
T0(60/dx+2:N) = kapa2 * dudx(60/dx+2:N);
dt = 0.01;

kFT = 1i * kx .* fft(T0);                    
v1 = v0 + dt * rho * real(ifft(kFT,'symmetric'));   
kFv = 1i * kx .* fft(v0);
dvdx = ifft(kFv, 'symmetric');
T1(1:60/dx+1) = T0(1:60/dx+1) + dt * kapa1 * real(dvdx(1:60/dx+1));
T1(60/dx+2:N) = T0(60/dx+2:N) + dt * kapa2 * real(dvdx(60/dx+2:N));

figure('Units','centimeter','Position',[2 2 25 40]);
subplot(11,2,1)
plot(lc,v0)
ylim([-0.3, 0.3])
title('velocity at time = 0 s');
subplot(11,2,2)
plot(lc,T0)
ylim([-0.3, 0.3])
title('stress at time = 0 s');

for i = 1:10000
    kFT = 1i * kx .* fft(T1);                    % ik F(T)
    v2 = v0 + 2 * dt * rho * real(ifft(kFT,'symmetric'));   % dv/dt = rho * dT/dx
    kFv = 1i * kx .* fft(v1);
    dvdx = ifft(kFv, 'symmetric');
    T2(1:60/dx+1) = T0(1:60/dx+1) + 2 * dt * kapa1 * real(dvdx(1:60/dx+1));
    T2(60/dx+2:N) = T0(60/dx+2:N) + 2 * dt * kapa2 * real(dvdx(60/dx+2:N));
    if mod(i,1000) == 0
        subplot(11,2,2*(i/1000)+1)
        plot(lc, real(v2))
        ylim([-0.3, 0.3])
        title(sprintf('velocity at time = %d s', i*dt))
        subplot(11,2,2*(i/1000)+2)
        plot(lc, real(T2))
        ylim([-0.3, 0.3])
        title(sprintf('stress at time = %d s', i*dt))
    end
    v0 = v1; v1 = v2;
    T0 = T1; T1 = T2;
end
xlabel('location(m)')
subplot(11,2,21)
xlabel('location(m)')
saveas(gcf, 'inhomo.png');

%% boundary conditions (Dirichlet conditions)
kapa = 1; rho = 1;
kx = fftshift((-N:N-1) * (2*pi / (2*N*dx))); 
v0 = zeros(1,2*N);
T0 = zeros(1,2*N);
T0(1:N) = FT; T0(N+1:2*N) = -FT;
dt = 0.01;

kFT = 1i * kx .* fft(T0);                    % ik F(T)
v1 = v0 + dt * rho * ifft(kFT,'symmetric');   % dv/dt = rho * dT/dx
kFv = 1i * kx .* fft(v0);
T1 = T0 + dt * kapa * ifft(kFv, 'symmetric');

figure('Units','centimeter','Position',[2 2 25 40]);
subplot(11,2,1)
plot(lc,v0(1:N))
ylim([-0.3, 0.3])
title('velocity at time = 0 s');
subplot(11,2,2)
plot(lc,T0(1:N))
ylim([-0.3, 0.3])
title('stress at time = 0 s');

for i = 2:10000
    kFT = 1i * kx .* fft(T1);                    % ik F(T)
    v2 = v0 + 2 * dt * rho * ifft(kFT,'symmetric');   % dv/dt = rho * dT/dx
    kFv = 1i * kx .* fft(v1);
    T2 = T0 + 2 * dt * kapa * ifft(kFv, 'symmetric');
    if mod(i,1000) == 0
        subplot(11,2,2*(i/1000)+1)
        plot(lc, real(v2(1:N)))
        ylim([-0.3, 0.3])
        title(sprintf('velocity at time = %d s', i*dt));
        subplot(11,2,2*(i/1000)+2)
        plot(lc, real(T2(1:N)))
        ylim([-0.3, 0.3])
        title(sprintf('stress at time = %d s', i*dt));
    end
    v0 = v1; v1 = v2;
    T0 = T1; T1 = T2;
end
xlabel('location(m)')
subplot(11,2,21)
xlabel('location(m)')
saveas(gcf, 'dirichlet.png');

%% boundary conditions (Neumann conditions)
kapa = 1; rho = 1;
kx = fftshift((-N:N-1) * (2*pi / (2*N*dx))); 
v0 = zeros(1,2*N);
T0 = zeros(1,2*N);
T0(1:N) = FT; T0(N+1:2*N) = FT;
dt = 0.01;

kFT = 1i * kx .* fft(T0);                    % ik F(T)
v1 = v0 + dt * rho * ifft(kFT,'symmetric');   % dv/dt = rho * dT/dx
kFv = 1i * kx .* fft(v0);
T1 = T0 + dt * kapa * ifft(kFv, 'symmetric');

figure('Units','centimeter','Position',[2 2 25 40]);
subplot(11,2,1)
plot(lc,v0(1:N))
ylim([-0.3, 0.3])
title('velocity at time = 0 s');
subplot(11,2,2)
plot(lc,T0(1:N))
ylim([-0.3, 0.3])
title('stress at time = 0 s');

for i = 2:10000
    kFT = 1i * kx .* fft(T1);                    % ik F(T)
    v2 = v0 + 2 * dt * rho * ifft(kFT,'symmetric');   % dv/dt = rho * dT/dx
    kFv = 1i * kx .* fft(v1);
    T2 = T0 + 2 * dt * kapa * ifft(kFv, 'symmetric');
    if mod(i,1000) == 0
        subplot(11,2,2*(i/1000)+1)
        plot(lc, real(v2(1:N)))
        ylim([-0.3, 0.3])
        title(sprintf('velocity at time = %d s', i*dt));
        subplot(11,2,2*(i/1000)+2)
        plot(lc, real(T2(1:N)))
        ylim([-0.3, 0.3])
        title(sprintf('stress at time = %d s', i*dt));
    end
    v0 = v1; v1 = v2;
    T0 = T1; T1 = T2;
end
xlabel('location(m)')
subplot(11,2,21)
xlabel('location(m)')
saveas(gcf, 'neumann.png');
