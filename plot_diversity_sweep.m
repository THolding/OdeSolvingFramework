name = "diversity_sweep_test";

numStrainsList = [];
prevalenceList = [];
prevalenceList2 = [];

runCount = 0;
while (exist([name num2str(runCount) '.csv'], 'file') == 2)
	numStrains = csvread([name num2str(runCount) "_numStrains.csv"]);
	numStrainsList = [numStrainsList, numStrains];
	data = csvread([name num2str(runCount) ".csv"]);

	t_index = 1;
	zs_start_index = t_index+1;
	zs_end_index = zs_start_index+numStrains-1;
	ws_start_index = zs_end_index+1;
	ws_end_index = ws_start_index+numStrains-1;
	ys_start_index = ws_end_index+1;
	ys_end_index = ys_start_index+numStrains-1;

	t = data(t_index,:);
	zs = data(zs_start_index:zs_end_index, :);
	ws = data(ws_start_index:ws_end_index, :);
	ys = data(ys_start_index:ys_end_index, :);

	graphics_toolkit('gnuplot')

	subplot(1,2,1);
	hold on;
	plot(t, zs, 'LineWidth', 2);
	%for i=1:(size(zs,1))
	%	plot(t, zs(i,:), 'b', 'LineWidth', 2);
	%end
	legend;
	xlabel('t');
	ylabel('z');

	subplot(1,2,2);
	hold on;
	plot(t, ys, 'LineWidth', 2);
	legend;
	xlabel('t');
	ylabel('y');

	hold off;
	print('-dpng', [name num2str(runCount) ".png"])
	close;

	%calculate prevalence
	
	prevalence = sum(ys(:,length(ys)));
	prevalenceList = [prevalenceList, prevalence];

	%Moving windowed mean:
	windowSize = 50;
	threshold = t(length(t)-windowSize);
	toMean = t>(t(length(t))-50);
	prevalence2 = sum(mean(ys(:,toMean),2));
	prevalenceList2 = [prevalenceList2, prevalence2];

	runCount = runCount+1;
end

figure
plot(numStrainsList, prevalenceList, 'k', 'LineWidth', 2);
hold on;
plot(numStrainsList, prevalenceList2, 'r', 'LineWidth', 2);
xlabel("number of strains");
ylabel("prevalence");

print('dpdf', [name "_prevalence.pdf"]);
