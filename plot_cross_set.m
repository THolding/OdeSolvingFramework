name = "si_anon_strains_cross";
numRuns = 5;

numStrainsList = csvread('si_anon_strains_numstrains_list.csv');

colours = [];
for (index = 1:numRuns)
	colours = [colours; index/numRuns, 0.4, 0];
end

hold on;
for (index = 1:numRuns)
	curName = [name, num2str(index-1), "_numstrains_sweep_prevalence.csv"];
	disp(curName);
	curPrevalence = csvread(curName);
	plot(numStrainsList, curPrevalence, 'LineWidth', 2, 'Color', colours(index,:));
end

name = "si_anon_strains_exp_cross";
colours = [];
for (index = 1:numRuns)
	colours = [colours; 0, 0.4, index/numRuns];
end

for (index = 1:numRuns)
	curName = [name, num2str(index-1), "_numstrains_sweep_prevalence.csv"];
	disp(curName);
	curPrevalence = csvread(curName);
	plot(numStrainsList, curPrevalence, 'LineWidth', 2, 'Color', colours(index,:));
end

legend(["0.00 (lin)"; "0.25 (lin)"; "0.50 (lin)"; "0.75 (lin)"; "1.00 (lin)"; "0.00 (exp)"; "0.25 (exp)"; "0.50 (exp)"; "0.75 (exp)"; "1.00 (exp)"], 'Location', 'SouthEast');
