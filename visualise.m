%        Author:  Ben Stott, bgbnbigben@contextualsystems.com
%  Organization:  Contextual Systems
%       Created:  26/12/12 21:48:40
% Last Modified:  Sat 29 Dec 2012 16:00:46 EST
%
%
% VISUALISE 
%   visualise(x, y, size, plotfunction) creates a Nelder-Mead visualisation,
%   drawn over a given plot function
%       It operates on the following input:
%           x a vector of x values!
%           y is a vector of y values
%           size is the size of the viewing window
%           plotfunction is a string to a matlab function to allow the x/y
%               values to be drawn over a coherent picture
%       
%
%       It has no functional output (as it draws to screen)
%       

function visualise(x, y, func, plotfunction, params)

    persistent initialised;
    persistent Frames; 
    persistent minval;
    persistent iter_no;
    
    win = params(1);
    vw = params(2:3);
    reruns = params(4);
    fps = params(5);
    max_steps = params(6);
    steps = 0:max_steps;
    
    if isempty(initialised)
        initialised = 1;
        title('Nelder-Mead animation')
        Frames = moviein(max_steps + 1);
        minval = zeros(1, max_steps + 1);
        iter_no = 1;
    else
        iter_no = iter_no + 1;
    end
    
    minval(iter_no) = feval(func, x( :, 1));
    A = [x(1, :), x(1, 1)]';
    B = [x(2, :), x(2, 1)]';
    C = [y, y(1)]';
    [X, Y] = meshgrid(-win:win/20:win);
    Z = feval(plotfunction, X, Y);
    subplot(1, 4, 1:3)
    mesh(X, Y, Z), xlabel('x'), ylabel('y'), zlabel('z'), view(vw);
    hold on;
    plot3(A, B, C, 'k', 'LineWidth', 2), title('Nelder-Mead Simplex Visualization');
    hold off;
    subplot(1, 4, 4)
    plot(steps, minval), title('Objective Function Value');
    Frames(:, iter_no) = getframe;

end
