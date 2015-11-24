data = csvread("repertoire_multistrain_w.csv");
numStrains = csvread("repertoire_multistrain_w_numStrains.csv");

t_index = 1;
zs_start_index = t_index+1;
zs_end_index = zs_start_index+numStrains-1;
ws_start_index = zs_end_index+1;
ws_end_index = ws_start_index+numStrains-1;
ys_start_index = ws_end_index+1;
ys_end_index = ys_start_index+numStrains-1;



t = data(t_index,:);
zs = data(zs_start_index:zs_end_index, :);
ws = data(ws_start_index:ws_end_index, :);
ys = data(ys_start_index:ys_end_index, :);

subplot(1,2,1);
plot(t, zs(1,:), 'b', 'LineWidth', 2);
hold on;
plot(t, zs(2,:), 'r', 'LineWidth', 2);
plot(t, zs(3,:), 'g', 'LineWidth', 2);
plot(t, zs(4,:), 'g', 'LineWidth', 2);
plot(t, zs(5,:), 'r', 'LineWidth', 2);
plot(t, zs(6,:), 'b', 'LineWidth', 2);
legend;
xlabel('t');
ylabel('z');

graphics_toolkit('gnuplot')

subplot(1,2,2);
plot(t, ys, 'LineWidth', 2);
legend;
xlabel('t');
ylabel('y');
