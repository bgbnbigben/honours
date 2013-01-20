%        Author:  Ben Stott, bgbnbigben@contextualsystems.com
%  Organization:  Contextual Systems
%       Created:  12/01/13 19:49:01
% Last Modified:  Sat 12 Jan 2013 21:36:01 EST
%
%
% ROSENBROCK 
%   [fx] = rosenbrock(x) computes .
%       It operates on the following input:
%       
%
%       The resultant output is:
%       

function [fx] = rosenbrock(x)
    x0 = x(1:end-1);
    x1 = x(2:end);

    fx = sum( (100*(x1-x0.^2).^2 + (x0-1).^2) , 2);
    %fx = (1 - x(1)).^2 + 100*(x(2) - x(1).^2).^2;
    %fx = x'*x;
end

