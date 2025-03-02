clc; clear; close all;

steps = load("list.txt");
n = length(steps);
dt = 0.6;
%% forward plotting

filename = "T_forward.h5";
figure('Position', [100, 100, 400, 1200]);

for i = 1:n
    time = num2str(steps(i)*dt, '%.2f');
    dataset_name = sprintf("/Time_%s", time);
    a = h5read(filename, dataset_name);
    
    subplot(n, 1, i);
    plot(0:100, a, 'k');
    text(5, max(a) * 0.9, sprintf('time = %.2f', steps(i) * dt), 'FontSize', 10, 'FontWeight', 'bold');
    if i < n
        set(gca, 'XTickLabel', []);
    end
    set(gca, 'Box', 'off');
end
filetitle = sprintf('Explicit scheme(dt = %.2fs)', dt);
sgtitle(filetitle);
filename = sprintf('forward_%.2f.png', dt);
saveas(gcf, filename);

%% backward plotting
filename = "T_backward.h5";
figure('Position', [100, 100, 400, 1200]);

for i = 1:n
    time = num2str(steps(i)*dt, '%.2f');
    dataset_name = sprintf("/Time_%s", time);
    a = h5read(filename, dataset_name);
    
    subplot(n, 1, i);
    plot(0:100, a, 'k');
    text(5, max(a) * 0.9, sprintf('time = %.2f', steps(i) * dt), 'FontSize', 10, 'FontWeight', 'bold');
    if i < n
        set(gca, 'XTickLabel', []);
    end
    set(gca, 'Box', 'off');
end
filetitle = sprintf('Implicit scheme(dt = %.2fs)', dt);
sgtitle(filetitle);
filename = sprintf('backward_%.2f.png', dt);
saveas(gcf, filename);

%% Crank-Nicolson Scheme
filename = "T_CNscheme.h5";
figure('Position', [100, 100, 400, 1200]);

for i = 1:n
    time = num2str(steps(i)*dt, '%.2f');
    dataset_name = sprintf("/Time_%s", time);
    a = h5read(filename, dataset_name);
    
    subplot(n, 1, i);
    plot(0:100, a, 'k');
    text(5, max(a) * 0.9, sprintf('time = %.2f', steps(i) * dt), 'FontSize', 10, 'FontWeight', 'bold');
    if i < n
        set(gca, 'XTickLabel', []);
    end
    set(gca, 'Box', 'off');
end
filetitle = sprintf('Semi-implicit scheme(dt = %.2fs)', dt);
sgtitle(filetitle);
filename = sprintf('CNscheme_%.2f.png', dt);
saveas(gcf, filename);
