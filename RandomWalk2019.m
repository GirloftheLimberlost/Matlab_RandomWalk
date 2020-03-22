% This code considers the random-walk motion of a single particle over the
% course of steps.
%if you don't see steps # of points @ the end on the final plot, then that
%means that the particle backtraced into one of the old spots
%preallocate global variables
counter = 1; %keeps track of steps taken
posnx1 = 0; %keeps track of coordinate x
posny1 = 0; %keeps track of cooordinate y
% steps vector defines different random walk trajectories
% used to create a 2-level nested for loop, which creates single plots for
% each
steps=[100,1000,5000,10000];
% HALFst(i) = 1/2 .* steps(i);
% will be used in the 2-level for loop;
HALFst= zeros(1,length(steps));
for i=1:length(steps);
HALFst(i)=steps(i)/2;
end
threshx = 0.5;
threshy = 0.5;
for i=1:length(steps)
for j = 1:steps(i)
% preallocate data_structure to store the 1- 10 trials of random walk
% for each N_step value
Distances_X = zeros (10,steps(i)); %holds the x distances
Distances_Y= zeros(10,steps(i)); % holds the y distances
d = zeros(10,steps(i)); %combines x and y distance
%plot(Distances_X, Distances_Y);
% gca = get current axes or figure
% set (gca, 'XLim', [lower_limit, upper_limit])
% the following two lines manually set the axes limits
% upper lim and lower lim = +/-(steps/2)
12
set(gca,'XLim',[-HALFst(i), HALFst(i)]);
set(gca,'YLim',[-HALFst(i), HALFst(i)]);
% for each step, get the corresponding frame
%F(j) = getframe;
for k=1:10
for s = 1:steps(i)
% temporary holder variable is a random number
% rand (rows, columns) [=] tmp has dimensions 2 rows, 1 column
tmp = rand(2,1);
if tmp(1,1) > threshx
posnx1 = posnx1 + 1;
else
posnx1 = posnx1 - 1;
end
if tmp(2,1) > threshy
posny1 = posny1 + 1;
else
posny1 = posny1 - 1;
end
% store tmp indices into Distances_X(i,k) and Distances_Y(i,k)
Distances_X(k,s) = posnx1;
Distances_Y(k,s) = posny1;
end
end
%Average values in the d vector which holds all the distances
%d_AVG= zeros(1,steps(i));
d = (Distances_X.^2 + Distances_Y.^2).^0.5;
d_AVG= mean(d); % mean (d) gives a row vector of the mean of each colummn
log_AVG = log(d_AVG);
log_steps = log(1:steps(i));
end
figure(i);
hold on;
%loglog(log_steps, log_AVG,'-s');
xlabel('log (N_{steps})');
ylabel('(Average displacement)');
title('Log-Log scale of Average Displacement vs. N = 10,000 steps');
legend('Random Walk points', 'Best fit estimation');
N=log_steps; % "x" value [=] log steps taken
c= log_AVG; % "y" value [=] log average distance given steps taken
p = polyfit(N,c,2);
%f(N) = cN^p
%G = log(f(N))= log(cN^p)
%G = log(c) + (p(3) .* log(N));

G=(polyval(p,c));

%newY=polyval(coeff2, log_x_value);

hold n; 
loglog(N,c, '*')

holdon; 
loglog(N,G, 'o');

hold off;

end
