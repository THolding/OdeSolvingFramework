numReps = 40;

graphics_toolkit('gnuplot')

for i=0:(numReps-1)
	name = ["sis_eirvsprev_" num2str(i)]
	data = csvread([name, ".csv"]);
	
	t = data(1,:);
	S = data(2,:);
	I = data(3,:);
	R = data(4,:);
	
	plot(t, S, 'g', 'LineWidth', 2);
	hold on;
	plot(t, I, 'r', 'LineWidth', 2);
	plot(t, R, 'b', 'LineWidth', 2);

	xlabel('t');
	ylabel('prevalence');
	hold off;
	close;

	print(["plots_check/plot_" num2str(i) ".png"], '-dpng'); 
end
