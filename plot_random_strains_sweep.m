name = "random_num_strains_sweep";

selectedStrainsList = csvread([name "_selectedStrains.csv"]);
prevalenceList = csvread([name "_prevalence.csv"]);


figure
plot(selectedStrainsList, prevalenceList, 'k', 'LineWidth', 2);
xlabel("# randomly selected strains");
ylabel("prevalence");

print('dpdf', [name "_prevalence.pdf"]);
