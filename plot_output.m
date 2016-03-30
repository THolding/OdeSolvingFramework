close all

data = csvread("sir_example.csv");

t = data(1,:);
S = data(2,:);
I = data(3,:);
R = data(4,:);

graphics_toolkit('gnuplot')

plot(t, S, 'g', 'LineWidth', 2);
hold on;
plot(t, I, 'r', 'LineWidth', 2);
plot(t, R, 'b', 'LineWidth', 2);

xlabel('t');
ylabel('proportion');
