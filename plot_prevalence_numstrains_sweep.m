%Use to plot output from simulation_sets::eir_sweep.
%function [] = plot_prevalence(datafile, indicesToPlot)
	datafile = "num_strains_sweep.csv";

	data = csvread(datafile);

	numStrains = data(:,1);
	prevalence = data(:,2);

	graphics_toolkit('gnuplot')

	plot(numStrains, prevalence, 'k', 'LineWidth', 3);

	xlabel('No. strains');
	ylabel('Prevalence');

	%axis([0, 750, 0, 1.0])

	%printing
	%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3])
	%print -dpdf "SIR_prev_fit.pdf" -r100

%end
