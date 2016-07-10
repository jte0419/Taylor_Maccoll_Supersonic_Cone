function [zdot] = TM_EQUATIONS(theta,z,gam)

% Taylor-Maccoll Equation in State-Space Form
% Written by: JoshTheEngineer
% YouTube   : www.youtube.com/JoshTheEngineer
% Website   : www.JoshTheEngineer.com
% Started: 01/14/16
% Updated: 01/14/16 - Started code
%                   - Works as intended
%          01/15/16 - Added comments
%
% PURPOSE
% - State-space form of the Taylor-Maccoll equation used for integration
%
% INPUTS
% - theta : Integration angle [rad]
% - z     : Angular and radial velocities []
% - gam   : Ratio of specific heats []
% 
% OUTPUTS
% - zdot  : Solution array

% Initialize solution vector
zdot = zeros(2,1);                                                          % zdot(1) = dVr/dTheta
                                                                            % zdot(2) = d2Vr/dTheta2
% Define term used often in equation below
A = (gam-1)/2;

% Numerator and denominator for zdot calculation below
num = (-2*A*z(1)) - (A*z(2)*cot(theta)) + (2*A*z(1)^3) + ...                % Numerator of second state-space equation
      (A*z(1)^2*z(2)*cot(theta)) + (2*A*z(1)*z(2)^2) + ...
      (A*z(2)^3*cot(theta)) + (z(1)*z(2)^2);
den = A*(1-z(1)^2-z(2)^2) - z(2)^2;                                         % Denominator of second state-space equation

% State-space representation
zdot(1) = z(2);
zdot(2) = num/den;