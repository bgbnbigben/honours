function v = simplex(type,x,b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMPLEX generates a matrix whose columns are the vertices
% of an orthogonal simplex with root-vertex x.
%
% v = simplex(type, x, b)
%
% Input type, x the roor vertex and b a scaling factor.
%   Returns a simplex vertices with maximal side length of 1
% Where type is either:
%       mat = Matlab's method based on initial vertex x
%       qr = QR decompostion of basis b about initial vertex x
%      ijk = orth. reg. simplex about x 
%       hh = orth. simplex about x which is the orth. complement of -b
%       random = uses rand function to generate simplex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2,
  error('Insufficient information supplied')
end
        
% initialise
n = size(x,1);
        
if n ==1
    x = x(:);
    n = length(x);
end
        
v = zeros(n,n+1);
v(:,1) = x(:,1);
        
        switch type
            case 'mat'
                usual_delta = 0.05;
                zero_term_delta = 0.00025;
                for j = 1:n
                    y = x(:,1);
                    if y(j) ~= 0
                        y(j) = (1 + usual_delta)*y(j);
                    else
                        y(j) = zero_term_delta;
                    end
                    v(:,j+1) = y;
                end
                
            case 'ijk'
                % create a regular orthogonal simplex
                b = b*ones(1,n);
                % create an orthogonal simplex using side_lengths
                % check there are no zero lengths
                for j = 1:n
                        y = x(:,1);
                        y(j) = y(j) + b(j);
                        v(:,j+1) = y;
                 end
     
            case 'hh'
                % create a simplex about x(:,1) with orthogonal
                % decomposition from vector -b
                % sum of orthogonal decompostion vectors = -b
                
                % given a simplex, find the longest basis vector
                [basis basis_lengths] = frame('pb',x);
                j = find(basis_lengths(1:n) == max(basis_lengths(1:n)));
                j = j(1);
                bl = basis(:,j);
                
                % calculate the orthogonal decompostion, vectors sum to -b
                basis_orth = hh(bl);
                
                % create the new simplex
                for j = 1:n
                    v(:,j+1) = v(:,1) + basis_orth(:,j);
                end
                
            case 'qr'
                %create simplex about x(:,1) using QR decompostion
                % of basis vectors for simplex x
                
                % get the basis vectors for the current simplex
                [basis basis_lengths] = frame('pb',x);
                
                % order basis vectors according to length of first
                % n basis vectors
                [sorted_lengths, j] = sort(basis_lengths(1:n));
                
                %get in descending order
                j = fliplr(j);
                basis = basis(:,j);
                
                % find QR decomposition of the ordered basis vectors
                [Q R] = qr(basis);
                
                % setup new length criteria
                d = diag(R);
                davg = sum(abs(d)) / n;
                sign_d = sign(0.5 +sign(d));
                d_new = sign_d .* max(abs(d), davg/10);
                D = diag(d_new);
                
                % calculate new basis vectors
                basis = Q*D;
                
                % create new simplex about x(:,1)
                % new simplex is x(:,1) and x(:,1) + basis(:,j) for j=1..n
                for j = 1:n
                    v(:,j+1) = v(:,1) + basis(:,j);
                end
            case 'random'
                v(:,1)=x(:,1);
                v(:,2:n+1)=x(:,1)*ones(1,n)+b*eye(n,n);
            otherwise
                error('An unknown simplex type has been used here')
        end
        
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
end

function B = hh(b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For a given vector b, B = HH(b) returns a matrix  B
% whose columns are orthogonal and sum(B:,i)) = -b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = length(b);
t = norm(b)/sqrt(n)*sign(sign(b(1)+.5));
v = b + t;
B = t*(eye(n)-((2/(v'*v))*v)*v');
end
end

