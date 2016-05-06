%Use to plot output from simulation_sets::eir_sweep.
%function [] = plot_prevalence(datafile, indicesToPlot)
	datafile = "eir_sweep.csv";

	data = csvread(datafile);

	eir = data(:,1);
	prevalence = data(:,indicesToPlot);

	%convert eir into annual not daily
	eir = eir*365;

	graphics_toolkit('gnuplot')

	plot(eir, prevalence, 'k', 'LineWidth', 3);

	xlabel('Annual EIR');
	ylabel('Prevalence');

	%axis([0, 750, 0, 1.0])

	%printing
	%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3])
	%print -dpdf "SIR_prev_fit.pdf" -r100

%end
