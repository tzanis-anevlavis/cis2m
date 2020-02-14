% Value function for the minimum time control of the double integrator.
% The set { x | V(x) <= T } is equal to the set of all states that can be
% steered to the origin in T seconds
function [ V ] = doubleIntegValueFun( x )
a = x(1) + 0.5*x(2)*abs(x(2));
if(a > 0)
    V = x(2) + 2*sqrt(x(1)+0.5*x(2)^2);
else
    V = -x(2) + 2*sqrt(-x(1)+0.5*x(2)^2);
end
end
