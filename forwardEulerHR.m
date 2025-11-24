function y = forwardEulerHR(func,t,y1)

% compute a forward Euler approximation to 
% dy/dt = f(t,y) with initial condition y(t1)=y1
% t is a vector of evenly spaced points: [t1,t2,...,tn]

% initialize y to be the same size at t
y = 0 * t;

% use the initial condition
y(1)=y1;

% use backward Euler to find y_i+1 given y_i
for i=1:length(y)-1
    % advance in time by evaluating y'(t) at the
    % current point ( t(i), y(i) )
    y(i+1) = y(i) + (t(i+1)-t(i))*feval(func,t(i),y(i));
end