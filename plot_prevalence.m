close all

data = csvread("eirvsprev_prevalences.csv");

eir = data(:,1);
prevalence = data(:,2);

%convert eir into annual not daily
eir = eir*365;

graphics_toolkit('gnuplot')

plot(eir, prevalence, 'g', 'LineWidth', 3);

xlabel('Annual EIR');
ylabel('prevalence');

axis([0, 750, 0, 1.0])

%printing
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3])
print -dpdf "SIRS_prev_fit.pdf" -r100

