function [diam basis Sdet norm_max] = Sdiam(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Returns the simplex diameter and basis vectors formed
%	from the sides of the simplex, the determinant and 
%   the maximum length of the basis vectors (simplex sides).
%   Input:   s the vertices of a simplex
%   Output:  diam - the diameter of the simplex
%            basis - the set of sides translating s_1 to the origin
%            Sdet - The determinant of the basis matrix
%            norm_max - the maximal length of the basis set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = length(s(:, 1));
basis = s(:, 2:n+1) - repmat(s(:, 1), 1, n);
Sdet = det(basis); 
d = 0;
norm_max = 0;
for k = 1:n
	b = norm(basis(:, k), 2);
	norm_max = max(norm_max, b);
	d = d + b; 
end
diam = d/n;
