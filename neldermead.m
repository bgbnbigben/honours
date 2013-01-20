function [xmin,ymin] = neldermead(n, f, xbar, rad, max_steps, g, params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: n - dimension of problem
%        f - function to be minimized
%        xbar - initial guess (column vector of length n) / centroid of simplex
%        rad - radius about initial guess
%        max_steps - maximum number of steps
%  
% Params:
%
%        g - plot version of function to be minimized
%
%
% Output: xmin - optimal solution
%         ymin - optimal objective value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize Simplex

x(:,1)=xbar;
x(:,2:n+1)=xbar*ones(1,n)+rad*eye(n,n);
for j=1:n+1
  y(j)=feval(f,x(:,j));
end  

% First frame of animation
visualise(x, y, f, g, params);

[y,r]=sort(y);
x=x(:,r);
diam =Sdiam(x);
tol=10^(-5);
i=0;
while i<=max_steps && diam >tol
   i=i+1;
   xbar = mean(x(:,1:n)')';
   xh = x(:,n+1);
   xr = 2*xbar - xh; yr = feval(f,xr);
   if yr < y(n)
      if yr < y(1)
         xe = 3*xbar - 2*xh; ye = feval(f,xe);
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
         xoc = 1.5*xbar - 0.5*xh; yoc = feval(f,xoc);
         if yoc < yr
            x(:,n+1) = xoc; y(n+1) = yoc;
         else  
            for j = 2:n+1
               x(:,j) = 0.5*x(:,1)+0.5*x(:,j); y(j) = feval(f,x(:,j));
            end   
         end
      else
         xic = 0.5*xbar + 0.5*xh; yic = feval(f,xic);
         if yic < y(n+1)
            x(:,n+1) = xic; y(n+1) = yic;
         else  
            for j = 2:n+1
               x(:,j) = 0.5*x(:,1)+0.5*x(:,j); y(j) = feval(f,x(:,j));
            end   
         end      
      end
   end
   [y,r] = sort(y);
   x = x(:,r);
   diam =Sdiam(x);
   minval(i) = feval(f,x(:,1));


% Current frame of animation
   pause(.05);
   visualise(x, y, f, g, params);
      
end  

xmin = x(:,1);
ymin = y(1);

%disp(' ');
%disp('Hit any key to PLAY the MOVIE'); 
%pause

% Now play the movie: 

%movie(Frames,reruns,fps);
