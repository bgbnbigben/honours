function z = griewangk(x)
    z = (x'*x)/4000.0 - prod(cos(x./sqrt((1:length(x))'))) + 1.0;
end
