%Thanks to Ian Grundy for these files.
% Plotting/Movie Parameters

win = 100;             % x, y plotted in the interval [-win,win]
view = [350,70];      % first = angle around from negative y, second = angle up from xy plane
reruns = 1;           % number of times the movie will run
framespersecond = 2;  % number of frames per second; small value = slower

% Nelder-Mead Method Parameters

dim = 2;                                                  % size of problem - number of independent variables                                          
xmin = zeros(dim,1);                                      % best vertex of final simplex
ymin = 0;                                                 % function value at best vertex
guess = [win/2+win/2*rand;win/2+win/2*rand];              % random initial guess for location of minimum
radius = 2.0;                                             % for generating the initial simplex about "guess" 
max_steps = 50;                                           % maximum number of steps


parameters = [win, view, reruns, framespersecond]

% Run the Nelder Mead Algorithm 

[xmin, ymin] = neldermead(dim, 'paraboloid', guess, radius, max_steps, 'paraboloidplot', parameters);

% Display the Solution

disp('best = ');
disp(xmin');
disp('value = ');
disp(ymin);
