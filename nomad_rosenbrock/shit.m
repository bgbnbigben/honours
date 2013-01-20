s = tf('s');
H = 1/(1 + s);
for i=1:1
    CL = H/(1+i*H);
    hold on;
    step(CL); % step response
    hold off;
end