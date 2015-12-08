numBetaDecaysList = csvread("sis_anon_strains_betadecays_list.csv");
numStrains = csvread("sis_anon_strains_numStrains.csv");

graphics_toolkit('gnuplot')

prevalenceList = [];
for i=0:(length(numBetaDecaysList)-1)
	data = csvread(["sis_anon_strains_betadecay_sweep_" num2str(i) ".csv"]);
	betaDecay = numBetaDecaysList(i+1);
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
	
	prevalence = Itot(length(Itot));
	prevalenceList = [prevalenceList, prevalence];
	
	plot(t, Ss, 'g', 'LineWidth', 1);
	hold on;
	plot(t, Is, 'r', 'LineWidth', 1);
	plot(t, Stot, 'g', 'LineWidth', 2);
	plot(t, Itot, 'r', 'LineWidth', 2);

	xlabel('t');
	ylabel('proportion');
	title(["betaDecay = " num2str(betaDecay) " prevalence = " num2str(prevalence)]);
	hold off;
end

figure
plot(numBetaDecaysList, prevalenceList, 'k', 'LineWidth', 3);
xlabel('transmission decay rate (alpha)');
ylabel('prevalence in population');
csvwrite('sis_anon_strains_betadecay_sweep_prevalence.csv', prevalenceList);
