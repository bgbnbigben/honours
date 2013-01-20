function [sol, eval, numfunc] = neldermead_positivebasis(f, xbar, tol, rad,...
max_steps, type, g, params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input: n - dimension of problem
%        f - function to be minimized
%        tol - accuracy of final point
%        xbar - initial guess (column vector of length n) / centroid of simplex
%        rad - radius about initial guess
%        max_steps - maximum number of steps
%
% Output: sol - solution vector
%         eval - function value at sol
%         numfunc - number of function evaluations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
n = length(xbar);
y = zeros(1, n + 1);
x(:, 1) = xbar;

% Initialize Simplex
if nargin < 5
    max_steps = 15000;
end
if nargin < 4
    rad = 0.1;
end
if nargin < 6
    type = 'ijk';
end
x = simplex(type, x, rad);
diam = Sdiam(x);

% God damnit matlab why can't you be normal for once?
% Using a java hashtable instead of a Map because of the dumb restrictions
% placed on parfor loops
cache = java.util.Hashtable;
%cache = containers.Map;

parfor i1 = 1 : n+1
    y(i1) = feval(f, x(:, i1));
    cache.put(mat2str(x(:, i1)), y(i1));
end
numfunc = n;
% First frame of animation
visualise(x, y, f, g, params);

% Sort x and y values by ascending y
[y, r] = sort(y);
x = x( :, r);
iter_no = 0;
eps_tol = min(diam^(1.5), tol);
good = true;

while (abs(y(n+1) - y(1)) > tol || diam > tol) && iter_no <= max_steps && good
    iter_no = iter_no +1;
    good = false;
    
    xbar = mean(x( :, 1:n), 2);
    xh = x( :, n+1);
    xr = 2*xbar - xh;
    if cache.containsKey(mat2str(xr))
        yr = cache.get(mat2str(xr))
    else
        yr = feval(f, xr); numfunc = numfunc + 1;
        cache.put(mat2str(xr), yr);
    end
    if yr < y(n)
        if yr < y(1)
            xe = 3*xbar - 2*xh;
            if cache.containsKey(mat2str(xe))
                ye = cache.get(mat2str(xe));
            else
                ye = feval(f, xe); numfunc = numfunc + 1;
                cache.put(mat2str(xe), ye);
            end
            if (ye < yr) && (ye < y(n+1) - eps_tol)
                x( :, n+1) = xe;
                y(n+1) = ye;
                good = true;
            elseif yr < y(n+1) - eps_tol
                x( :, n+1) = xr;
                y(n+1) = yr;
                good = true;
            end
        elseif yr < y(n+1) - eps_tol
            x( :, n+1) = xr;
            y(n+1) = yr;
            good = true;
        end
    else
        if yr < y(n+1) - eps_tol
            xoc = 1.5*xbar - 0.5*xh;
            if cache.containsKey(mat2str(xoc))
                yoc = cache.get(mat2str(xoc));
            else
                yoc = feval(f, xoc); numfunc = numfunc + 1;
                cache.put(mat2str(xoc), yoc);
            end
            good = true;
            if yoc < yr
                x( :, n+1) = xoc;
                y(n+1) = yoc;
            else
                firstEl = x( :, 1);
                parfor j1 = 2 : n+1
                    x( :, j1) = 0.5*firstEl + 0.5*x( :, j1);
                    if cache.containsKey(mat2str(x( :, j1)))
                        y(j1) = cache.get(mat2str(x( :, j1)))
                    else
                        y(j1) = feval(f, x( :, j1));
                        cache.put(mat2str(x(:, j1)), y(j1));
                    end
                end
                numfunc = numfunc + n - 1;
            end
        else
            xic = 0.5 * xbar + 0.5 * xh;
            if cache.containsKey(mat2str(xic))
                yic = cache.get(mat2str(xic))
            else
                yic = feval(f, xic); numfunc = numfunc + 1;
                cache.put(mat2str(xic), yic);
            end
            if yic < y(n+1) - eps_tol
                x( :, n+1) = xic;
                y(n+1) = yic;
                good = true;
            end
        end
    end
    
    if ~good % we form a minimal frame and look to see if we have a quai-minimal point
        [diam basis] = Sdiam(x);
        temp = sum(basis( :, 1:n));
        if norm(temp, 2) < tol || det(basis) < tol
            x1 = x( :, 1);
            x = simplex('random', x, diam);
            x( :, 1) = x1;
            x = simplex(type, x, diam);
        else
            x( :, n+2) = x( :, 1) - diam * (temp(:) / norm(temp, 2));
            if cache.containsKey(mat2str(x( :, n+2)))
                y(n+2) = cache.get(mat2str(x( :, n+2)))
            else
                y(n+2) = feval(f, x( :, n + 2)); numfunc = numfunc + 1;
                cache.put(mat2str(x( :, n + 2)), y(n+2));
            end
            if y(n+2) < y(n+1) - eps_tol
                x( :, n+1) = x( :, n+2);
                y(n+1) = y(n+2);
                good = true;
                x( :, n+2) = [];
                y(n+2) = [];
            else
                firstEl = x(:, 1);
                parfor j1 = 2 : n+1
                    x( :, j1) = 0.5 * firstEl + 0.5 * x( :, j1);
                    if cache.containsKey(mat2str(x( :, j1)))
                        y(j1) = cache.get(mat2str(x( :, j1)))
                    else
                        y(j1) = feval(f, x( :, j1)); numfunc = numfunc + 1;
                        cache.put(mat2str(x( :, j1)), y(j1));
                    end
                end
                good = true;
                x( :, n+2) = [];
                y(n+2) = [];
            end
        end
    end
    
    [y,r] = sort(y);
    x = x( :, r);
    diam = Sdiam(x);
    eps_tol = min(tol, diam^(1.5));

% Current frame of animation
    visualise(x, y, f, g, params);
 
end

sol = x( :, 1);
eval = feval(f, sol); numfunc = numfunc + 1;
toc
keys(cache)
for key in cache.keys()
    key
end
end
