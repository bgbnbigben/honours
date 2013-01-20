function [F,G] = frame(type,v,h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The output from frame is determined by the type of
% infomation required.
% Code developed by David Byatt in his 2000 thesis
% Convergent variants of the Nelder-Mead Algorithm
% University of Canterbury,  NZ.
%
%  v is a matrix whose columns represent the vertices of
%  the simplex frame, h is a scale factor and type is either:
%  l = returns the n side lengths for the simplex v, scaled by h
%  f = complete the frame for the current simplex where the length
%      of the new frame point is scaled by h
% pb = return the n+1 positive basis vectors and their lengths
% sh = shrink the current frame v towards v(:,1) by the scale factor h
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2,
    error('Incorrect number of input arguments'),
end

% dimension
n = size(v,1);

if nargin == 2,
    h = 1;
end

% initialise
vectors = zeros(n,n+2);
lenghts = zeros(1,n+2);
G = [];

switch type
    case 'f'
        % return the completed frame
        vectors(:,1:n+1) = v;
        vectors(:,n+1) = (1 + h)*v(:,1) - h/n * sum(v(:,2:n+1),2);
        F = vectors;
    
    case 'sh'
        if size(v,2) ~= n+1,
            error('The input frame has the wrong dimensions.');
        end
        % shrink the current frame towards v(:,1) by h
        vectors(:,1) = v(:,1);
        for j = 2:n+1,
            vectors(:,j) = v(:,1) + h*(v(:,j) - v(:,1));
        end
        if vectors == v
            % changes are beyond machine precision
            vectors(:,2:n+2) = vectors(:,ones(1,n+1));
        end
        F = vectors;
    otherwise
        for j = 1:n,
            vectors(:,j) = (v(:,j+1) - v(:,1)) / h;
            lenghts(:,j) = norm(vectors(:,j));
        end
        switch type
            case 'pb'
                vectors(:,n+1) = -sum(vectors(:,1:n), 2) / n;
                lengths(:,n+1) = norm(vectors(:,n+1));
                F = vectors(:,1:n+1);
                G = lengths(:,1:n+1);
            case 'l'
                F = lengths(1:n);
            otherwise
                error('An unknown frame argument has been used')
        end
end
        