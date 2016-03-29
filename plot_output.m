name = "diversity_sweep_test0";


numStrains = csvread([name "_numStrains.csv"]);
data = csvread([name ".csv"]);

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

graphics_toolkit('gnuplot')

subplot(1,2,1);
hold on;
plot(t, zs, 'LineWidth', 2);
%for i=1:(size(zs,1))
%	plot(t, zs(i,:), 'b', 'LineWidth', 2);
%end
legend;
xlabel('t');
ylabel('z');

subplot(1,2,2);
hold on;
plot(t, ys, 'LineWidth', 2);
legend;
xlabel('t');
ylabel('y');

hold off;
