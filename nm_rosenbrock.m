% Thanks to Ian Grundy for these files.
% Plotting/Movie Parameters

win = -1.5;           % x, y plotted in the interval [-win,win]
view = [350,60];      % first = angle around from negative y, second = angle up from xy plane
reruns = 1;           % number of times the movie will run
framespersecond = 2;  % number of frames per second; small value = slower

% Nelder-Mead Method Parameters

dim = 2;                                           % size of problem - number of independent variables
guess = [-win/4-win*rand/2;win/4+win*rand/2];      % random initial guess for location of minimum
radius = 0.25;                                     % for generating the initial simplex about "guess" 
max_steps = 10000;                                 % maximum number of steps

parameters = [win, view, reruns, framespersecond, max_steps];

% Run the Nelder Mead Algorithm 
%%

%[ best vertex of final simplex, function value at best vertex ]

%[xmin, ymin] = neldermead(dim, 'rosenbrock', guess, radius, max_steps, 'rosenbrockplot', parameters);
[xmin, ymin, numfunc] = neldermead_positivebasis('rosenbrock', guess, 1e-100, radius, max_steps, 'ijk', 'rosenbrockplot', parameters);

% Display the Solution

fprintf('best x = [%s]\n', strtrim(sprintf('%f ', xmin')));
fprintf('corresponding y = %f\n', ymin);
fprintf('%d iterations\n', numfunc);
