%
% Heuristic Methods - Nelder-Mead Search
%
% Input: n - dimension of problem
%        f - function to be minimized
%        xbar - initial guess (column vector of length n) / centroid of simplex
%        rad - radius about initial guess
%        max_steps - maximum number of steps
%  
% Output: sol - solution vector 
%         eval - function value at sol
%         numfunc - number of function evaluations
%

function [sol, eval, numfunc] = neldermead_basic_ian(f, xbar, tol, rad, max_steps, type) 

numfunc = 0;
n=length(xbar);

% Initialize Simplex
if nargin <6
    max_steps = 15000;
end
if nargin <5
    rad = 0.5;
end
if nargin<7
    type ='ijk';
else
    x=simplex('random',xbar,rad);
end
x(:,1)=xbar;
x=simplex(type,x,rad);
diam = Sdiam(x);
y(1:n+1)=feval(f, x(:,1:n+1)); numfunc=numfunc+n+1;

[y,r]=sort(y);
x=x(:,r);

for i=1:max_steps

   xbar = mean(x(:,1:n)')';
   xh = x(:,n+1);
   xr = 2*xbar - xh; yr = feval(f,xr); numfunc=numfunc+1;
   if yr < y(n)
      if yr < y(1)
         xe = 3*xbar - 2*xh; ye = feval(f, xe); numfunc=numfunc+1;
         if ye < yr
            x(:,n+1) = xe; y(n+1) = ye;
         else  
            x(:,n+1) = xr; y(n+1) = yr;
         end
      else
         x(:,n+1) = xr; y(n+1) = yr;
      end
   else         
      if yr < y(n+1)
         xoc = 1.5*xbar - 0.5*xh; yoc = feval(f, xoc); numfunc=numfunc+1;
         if yoc < yr
            x(:,n+1) = xoc; y(n+1) = yoc;
         else  
            [x y numfunc] = shrink(f, x, numfunc, n);  
         end
      else
         xic = 0.5*xbar + 0.5*xh; yic = feval(f, xic); numfunc=numfunc+1;
         if yic < y(n+1)
            x(:,n+1) = xic; y(n+1) = yic;
         else  
            for j1 = 2:(n+1)
                x(:,j1) = 0.5*x(:,1)+0.5*x(:,j1); y(j1) = feval(f, x(:,j1)); numfunc=numfunc+1;
            end  
         end      
      end
   end
   [y,r] = sort(y);
   x = x(:,r);
         
end  

sol = x(:,1);
eval = feval(f, sol); 
numfunc=numfunc+1;

function [x y numfunc] = shrink(f, x, numfunc, n)
    for j1 = 2:(n+1)
        x(:,j1) = 0.5*x(:,1)+0.5*x(:,j1); y(j1) = feval(f, x(:,j1)); numfunc=numfunc+1;
    end
end
end


