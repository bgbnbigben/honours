function z = griewangk(x)
    %x = round(x);
    % Ackley's function
    %z = -20*exp(-0.2 *((1/10*(x'*x))^.5)) - exp(1/10 * (cos(2*pi*x(1)) + cos(2*pi*x(2))+ cos(2*pi*x(3))+ cos(2*pi*x(4))+ cos(2*pi*x(5))+ cos(2*pi*x(6))+ cos(2*pi*x(7))+ cos(2*pi*x(8))+ cos(2*pi*x(9))+ cos(2*pi*x(10)))) + 20 + exp(1);
    
    z = (x'*x)/4000.0 - prod(cos(x./sqrt((1:length(x))'))) + 1.0;
end
