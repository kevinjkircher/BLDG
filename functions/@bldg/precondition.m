function x0 = precondition(b,t,weather,gains,N,Ts)
%PRECONDITION finds a plausible initial state for the simplest building by
%simulating the building under the previous several days of weather and
%gains, and assumping perfect control.

% Syntax:
%	x0 = precondition(b,t,weather,gains,N,Ts)
%
% Inputs:
%	b, a bldg object containing building data.
%	t, the simulation time span in solar seconds since midnight on Jan. 1.
%	weather, a struct containing weather data.
%	gains, a struct containing internal gains data.
%	N, the number of finite difference nodes in the wall.
%	Ts, the indoor temperature setpoint.
%
% Outputs:
%	x0, the initial state of the building.

%
% Set the number of preconditioning days.
%
  ndShift = 14;
%
% Decide whether to use the previous or next days.
%
  dt = t(2) - t(1);
  if t(1) - ndShift*24*3600 > 0
	t0 = t(1) - ndShift*24*3600;			% use previous days if possible
    tf = t(1);
  else
	t0 = t(1);								% use next days if necessary
    tf = t(1) + ndShift*24*3600;
  end
  tp = t0:dt:tf;
%
% Extract the disturbances for preconditioning.
%
  Wp = getDisturbances(weather,gains,tp);
%
% Simulate the building. Start it at the temperature setpoint. Assume
% perfect control.
%
  xp0 = Ts*ones(N+2,1);
  [Xp,~] = bsim(b,tp,Wp,xp0);
  x0 = Xp(end,:)';
%
end

% Copyright 2016 Kevin J. Kircher. See bldg/license.txt for full license
% information.