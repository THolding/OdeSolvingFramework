name = 'numstrains_sweep';
%recovs = csvread([name, '_recovered.csv']);
%function [] = plot_eirprev_numstrains (name)

	numStrainsList = csvread([name, '_numStrainsList.csv']);
	data = csvread([name, '.csv'])';
	eir = data(1,:) * 365;
	numPrevs = length(numStrainsList);
	prevs = data(2:numPrevs+1, :);

	lineColours = calc_colour_range(1:length(numStrainsList), [.1,.2,.9], [0.9,0.2,0.1]);

	%Plot each eir sweep
	graphics_toolkit('gnuplot')
	hold off;
	for i=1:numPrevs
		plot(eir, prevs(i,:), 'Color', lineColours(i,:), 'LineWidth', 3);
		hold on;
	end

	legend(prepare_legend(numStrainsList, 'n='), 'Location', 'northwest');
	xlabel('EIR (annual)');
	%xlabel('beta');
	ylabel('Prevalence');
	%axis([0, eir(length(eir)), 0, ceil(max(prevs)*10)/10]); %Rounds axis scaling to 1/10th.

	%figure
	%plot(eir, recovs, 'b');

	%printing
	set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 4.5])
	print -dpdf "prevalence_vs_antigenic_diversity.pdf" -r100
%end