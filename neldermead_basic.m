function [sol, eval, numfunc] = neldermead_basic(f, xbar, tol, rad, max_steps, type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input: n - dimension of problem
%        f - function to be minimized
%        xbar - initial guess (column vector of length n) / centroid of simplex
%        rad - radius about initial guess
%        max_steps - maximum number of steps
%  
% Output: sol - solution vector 
%         eval - function value at sol
%         numfunc - number of function evaluations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
numfunc = 0;
n=length(xbar);

% Initialize Simplex
if nargin <5
    max_steps = 15000;
end
if nargin <4
    rad = 0.1;
end
if nargin<6
    type ='ijk';
else
    x=simplex('random',xbar,rad);
end
x(:,1)=xbar;
x=simplex(type,x,rad);
diam = Sdiam(x);
for i1=1:n+1
    y(i1)=feval(f, x(:,i1)); numfunc=numfunc+1;
end


[y,r]=sort(y);
x=x(:,r);
iter_no = 0;

while (abs(y(n+1) - y(1)) > tol ||diam > tol) && iter_no <= max_steps
iter_no = iter_no +1;

   xbar = mean(x(:,1:n),2);
   xh = x(:,n+1);
   xr = 2*xbar - xh; yr = feval(f, xr); numfunc=numfunc+1;
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
            for j1 = 2:n+1
                x(:,j1) = 0.5*x(:,1)+0.5*x(:,j1); y(j1) = feval(f, x(:,j1)); numfunc=numfunc+1;
             end   
         end
      else
         xic = 0.5*xbar + 0.5*xh; yic = feval(f, xic); numfunc=numfunc+1;
         if yic < y(n+1)
            x(:,n+1) = xic; y(n+1) = yic;
         else  
            for j1 = 2:n+1
                x(:,j1) = 0.5*x(:,1)+0.5*x(:,j1); y(j1) = feval(f, x(:,j1)); numfunc=numfunc+1;
            end 
         end      
      end
   end
   [y,r] = sort(y);
   x = x(:,r);
   diam = Sdiam(x);    
end  

sol = x(:,1);
eval = feval(f, sol); 
numfunc=numfunc+1;
toc
end
