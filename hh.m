function B = hh(b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For a given vector b, B = HH(b) returns a matrix  B
% whose columns are orthogonal and sum(B:,i)) = -b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = length(b);
t = norm(b)/sqrt(n)*sign(sign(b(1)+.5));
v = b + t;
B = t*(eye(n)-((2/(v'*v))*v)*v');