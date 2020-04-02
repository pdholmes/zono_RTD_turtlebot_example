function [] = generate_turtlebot_example_dynamics()
% dynamics -- x and y position, heading, velocity, acceleration, are the
% states, and acceleration and yaw rate (both assumed constant over this
% time horizon) are the trajectory parameters.

   syms x y theta v a w; % variables for FRS
   dx = v*cos(theta);
   dy = v*sin(theta);
   dtheta = w;
   dv = a;
   da = 0;
   dw = 0;
   
   z = [x; y; theta; v; a; w];
   dz = [dx; dy; dtheta; dv; da; dw];
   
   syms tdummy udummy % dummy variables
   
   matlabFunction(dz, 'File', 'dyn_turtlebot_example', 'vars', {tdummy z udummy});

end