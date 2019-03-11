function dy = lorenz(t,y,p1,p2)
% LORENZ   ODE for the Lorenz system
%     Usage with ode45:
%        [t x] = ode45('lorenz',[0 200],rand(1,3));
%     Changing default parameter r
%         r = 210;
%        [t x] = ode45('lorenz',[0 200],rand(1,3),[],r);

% parameter r 
if nargin<4
   r=28;
else
   r=p2;
end

dy = zeros(3,1);

% ODE of Lorenz 63
dy(1) = 10 * (y(2) - y(1));
dy(2) = r*y(1) - y(2) - y(1)*y(3);
dy(3) = y(1)*y(2) - (8/3)*y(3);

