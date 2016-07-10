function [sol] = THETA_BETA_M_v2(thetaInp,betaInp,MInp,gam,degrad)

% Written by: JoshTheEngineer
% YouTube   : www.youtube.com/JoshTheEngineer
% Website   : www.JoshTheEngineer.com
% Started: 01/14/16
% Updated: 01/14/16 - Started code
%                   - Works as intended
%          01/15/16 - Added some comments
%          01/17/16 - Incorporates different methods of solution
%                   - Allows for input/output of degrees or radians
% 
% PURPOSE
% - Use the Theta-Beta-M relation to obtain the flow deflection angle from
%   the shock wave angle and freestream Mach number
% - Set the input you wanted solved for to zero (0)
% - Make sure to specify whether the input angles are in 'deg' or 'rad'
% 
% REFERENCES
% - Modern Compressible Flow, Anderson, pg. 136, eqn. 4.17
%
% INPUTS
% - thetaInp : Cone angle to horizontal [deg/rad] 
% - betaInp  : Shockwave angle to horizontal [deg/rad]
% - MInp     : Freestream Mach number []
% - gam      : Ratio of specific heats, or cp/cv []
% - degrad   : Specify whether calculations are in degrees or radians [str]
%
% OUTPUTS
% - sol      : Solution value, depending on what is being solved for

% Convert between degrees and radians if necessary
if (strcmp(degrad,'deg'))
    thetaInp = thetaInp*(pi/180);
    betaInp  = betaInp*(pi/180);
elseif (strcmp(degrad,'rad'))
    % Do nothing, computations are in radians
end

% -------------------- SOLVING FOR: THETA ---------------------------------
if (thetaInp == 0)
    
    % Check against Mach wave angle equation
    if (betaInp <= asin(1/MInp))
        sol = 0;
        return;
    end
    
    % Numerator and denominator
    N = ((MInp*sin(betaInp))^2)-1;
    D = (MInp^2*(gam+cos(2*betaInp)))+2;
    
    % Wedge angle from the Theta-Beta-M relation
    theta = atan(2*cot(betaInp)*(N/D));
    
    % Set solution variable based on input deg/rad
    if (strcmp(degrad,'deg'))                                               % If the output is in degrees
        sol = theta*(180/pi);
    elseif (strcmp(degrad,'rad'))                                           % If the output is in radians
        sol = theta;
    end

% ----------------- SOLVING FOR: MACH NUMBER ------------------------------
elseif (MInp == 0)
    
    % Find zero of theta-beta-M equation
    myfun = @(t,b,g,M) (2*cot(b)*((((M*sin(b))^2)-1)/...
                        ((M^2*(g+cos(2*b)))+2))) - tan(t);
    t   = thetaInp;                                                         % Given cone angle [rad]
    b   = betaInp;                                                          % Given shock angle [rad]
    g   = gam;                                                              % Given ratio of specific heats []
    fun = @(M) myfun(t,b,g,M);                                              % Set given variables in myfun
    M   = fzero(fun,1);                                                     % Solve for M, starting at M = 1
    
    % Set solution variable
    sol = M;                                                                % Set the solution variable
    
% -------------------- SOLVING FOR: BETA ----------------------------------
elseif (betaInp == 0)
    
    % Find zero of theta-beta-M equation
    myfun = @(t,b,g,M) (2*cot(b)*((((M*sin(b))^2)-1)/...
                        ((M^2*(g+cos(2*b)))+2))) - tan(t);
    t   = thetaInp;                                                         % Given cone angle [rad]
    M   = MInp;                                                             % Given Mach number [rad]
    g   = gam;                                                              % Given ratio of specific heats []
    fun = @(b) myfun(t,b,g,M);                                              % Set given variables in myfun
    b   = fzero(fun,0.5);
    
    % Set solution variable based on input deg/rad
    if (strcmp(degrad,'deg'))                                               % If the output is in degrees
        sol = b*(180/pi);
    elseif (strcmp(degrad,'rad'))                                           % If the output is in radians
        sol  = b;
    end
    
end

