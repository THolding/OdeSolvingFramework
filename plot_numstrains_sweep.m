runName = "si_anon_strains_cross05_numstrains_sweep";
numStrainsList = csvread("si_anon_strains_numstrains_list.csv");

prevalenceList = [];
for i=0:(length(numStrainsList)-1)
	disp([runName "_" num2str(i) ".csv"])
	data = csvread([runName "_" num2str(i) ".csv"]);
	numStrains = numStrainsList(i+1);

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
	Stot = sum(Ss(1:numStrains, :));
	Rtot = Ss(numStrains+1, :); %== R
	Itot = sum(Is);
	
	prevalence = Itot(length(Itot));
	prevalenceList = [prevalenceList, prevalence];
	
	%plot(t, Ss, 'g', 'LineWidth', 1);
	%hold on;
	%plot(t, Is, 'r', 'LineWidth', 1);
	%plot(t, R, 'b', 'LineWidth', 1);
	%plot(t, Stot, 'g', 'LineWidth', 2);
	%plot(t, Itot, 'r', 'LineWidth', 2);
	%plot(t, Rtot, 'b', 'LineWidth', 2);

	%xlabel('t');
	%ylabel('proportion');
	%title(["n = " num2str(numStrains) " prevalence = " num2str(prevalence)]);
	%hold off;
end

figure
plot(numStrainsList, prevalenceList, 'k', 'LineWidth', 2);
xlabel('number of strains');
ylabel('prevalence in population');
csvwrite([runName '_prevalence.csv'], prevalenceList);


print('-dpdf', [runName '.pdf'])
