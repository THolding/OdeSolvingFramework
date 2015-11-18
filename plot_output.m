data = csvread("gupta_multistrain_w.csv");

t = data(1,:);
zs = data(2:5, :);
ws = data(6:9, :);
ys = data(10:13, :);

subplot(1,2,1);
plot(t, zs, 'LineWidth', 2);
legend;
xlabel('t');
ylabel('z');

subplot(1,2,2);
plot(t, ys, 'LineWidth', 2);
legend;
xlabel('t');
ylabel('y');
