data = csvread("sis_anon_strains.csv");
numStrains = csvread("sis_anon_strains_numStrains.csv");

t_index = 1;
Ss_start_index = t_index+1;
Ss_end_index = Ss_start_index+numStrains-1;
Is_start_index = Ss_end_index+1;
Is_end_index = Is_start_index+numStrains-1;

t = data(t_index,:);
Ss = data(Ss_start_index:Ss_end_index, :);
Is = data(Is_start_index:Is_end_index, :);
Stot = sum(Ss(1:numStrains, :));
Itot = sum(Is);

graphics_toolkit('gnuplot')

plot(t, Ss, 'g', 'LineWidth', 1);
hold on;
plot(t, Is, 'r', 'LineWidth', 1);
plot(t, Stot, 'g', 'LineWidth', 3);
hold on;
plot(t, Itot, 'r', 'LineWidth', 3);
xlabel('t');
ylabel('proportion');
