function [value,isterminal,direction] = EVENTS(theta,z,gam)

% Event Handling for ode15s Solver
% Written by: JoshTheEngineer
% YouTube   : www.youtube.com/JoshTheEngineer
% Website   : www.JoshTheEngineer.com
% Started: 01/14/16
% Updated: 01/14/16 - Started code
%                   - Works as intended
% 
% PURPOSE
% - Stops the ODE solution when the cone wall is reached
% - Only a single event is included right now, can expand as necessary
% 
% INPUTS
% - theta : Integration angle [rad]
% - z     : Angular and radial velocities []
% - gam   : Ratio of specific heats []
% 
% OUTPUTS
% - value      : Value of the ith event function
%                Event is triggered when value is 0
% - isterminal : Integration terminates at a zero of the event function
%                  if this value is 1, otherwise it's 0
% - direction  : If all zeros are to be located, 0
%                If zeros where event fnc decreasing, -1
%                If zeros where event fnc increasing, +1

value      = 1;             % No event triggered yet
isterminal = 1;             % Terminates after the first event is reached
direction  = 0;             % Get zeros from either direction (inc or dec)

% Trigger when angular velocity switches sign
% - Angular velocity needs to be zero at the cone surface (no suction or
%     blowing)
if (z(2) > 0)
    value = 0;              % If angular velocity z(2) crosses from
end                         %  being negative to positive, trigger event



