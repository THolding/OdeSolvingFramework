filen = 'numstrains_sweep_30_single_1.7500000000000009.csv';

disp(filen)
data = csvread(filen);

numStrains = (size(data,1)-2)/2;
t_index = 1;
s_start = t_index+1;
s_end = s_start + numStrains-1;
i_start = s_end+1;
i_end = i_start+numStrains-1;
r_index = i_end+1;

t = data(t_index,:);
Ss = data(s_start:s_end,:);
Is = data(i_start:i_end,:);
R = data(r_index,:);

%Calculate totals
sTot = sum(data(s_start:s_end,:),1);
iTot = sum(data(i_start:i_end,:),1);

%plot
figure
hold on;
plot(t, Ss, 'g', 'LineWidth', 1);
plot(t, Is, 'r', 'LineWidth', 1);
plot(t, R, 'b', 'LineWidth', 3);
plot(t, sTot, 'g', 'LineWidth', 3);
plot(t, iTot, 'i', 'LineWidth', 3);
xlabel('t');
ylabel('Proportion of population');
