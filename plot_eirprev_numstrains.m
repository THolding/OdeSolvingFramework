%Use to plot output from simulation_sets::num_strains_sweep_each_beta0.

name = 'num_strains_sweep_each_beta0';
%function [] = plot_eirprev_numstrains (name)

	beta0List = csvread([name, '_beta0List.csv']);
	numPrevs = length(beta0List);

	data = csvread([name, '.csv'])';
	numStrains = data(1,:);
	prevs = data(2:numPrevs+1, :);

	lineColours = calc_colour_range(1:length(beta0List), [.1,.2,.9], [0.9,0.2,0.1]);

	%Plot each prevalence vector against numStrains
	graphics_toolkit('gnuplot')
	hold off;
	for i=1:numPrevs
		plot(numStrains, prevs(i,:), 'Color', lineColours(i,:), 'LineWidth', 3);
		hold on;
	end

	legend(prepare_legend(beta0List, 'beta0 = '), 'Location', 'northwest');
	xlabel('No. strains');
	ylabel('Prevalence');
	%axis([0, eir(length(eir)), 0, ceil(max(prevs)*10)/10]); %Rounds axis scaling to 1/10th.

	%printing
	set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 4.5])
	print -dpdf "anon_sir_numstrains_vs_prevalence.pdf" -r100
%end
