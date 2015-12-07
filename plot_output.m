data = csvread("si_anon_strains.csv");
numStrains = csvread("si_anon_strains_numStrains.csv");

t_index = 1;
Ss_start_index = t_index+1;
Ss_end_index = Ss_start_index+numStrains;
Is_start_index = Ss_end_index+1;
Is_end_index = Is_start_index+numStrains-1;
R_index = Ss_end_index;



t = data(t_index,:);
Ss = data(Ss_start_index:Ss_end_index, :);
Is = data(Is_start_index:Is_end_index, :);
R = data(R_index, :);

graphics_toolkit('gnuplot')

plot(t, Ss, 'g', 'LineWidth', 2);
hold on;
plot(t, Is, 'r', 'LineWidth', 2);
plot(t, R, 'b', 'LineWidth', 2);

%legend;
xlabel('t');
ylabel('proportion');

figure;
Stot = sum(Ss(1:numStrains, :));
Rtot = Ss(numStrains+1, :); %== R
Itot = sum(Is);

plot(t, Stot, 'g', 'LineWidth', 2);
hold on;
plot(t, Itot, 'r', 'LineWidth', 2);
plot(t, Rtot, 'b', 'LineWidth', 2);
xlabel('t');
ylabel('proportion');
