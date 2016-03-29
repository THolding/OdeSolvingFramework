data = csvread('test.csv');

t = data(1,:);
Z = data(1,:);
X = data(1,:);
Y = data(1,:);

plot(t, Z, 'k');
hold on;
plot(t, X, 'b');
plot(t, Y, 'r');