function tspan = getTiming(tw,nd0,t0,tf)
%GETTIMING extracts the simulation time span from the weather measurement
%times.

% Inputs:
%	tw, an array containing a year of measurement times.
%	nd0, the initial simulation day number.
%	t0, the initial simulation time (seconds since midnight on day nd0).
%	tf, the final simulation time (seconds since midnight on day nd0).
%
% Output:
%	tspan, the simulation time span in seconds since midnight on day 1.
%
% The output time span has the same sample time as the measurements.

%
% Check that inputs make sense.
%
  if tf <= t0, error('Final time must exceed initial time.'), end
  if nd0+tf/24/3600 > 366
    error('Final time must be before midnight on December 31.')
  end
  if nd0+t0/24/3600 < 1
    error('Initial time must be after midnight on January 1.')
  end
%
% Get the sample time and number of output time steps.
%
  dt = tw(2) - tw(1);
  nt = floor((tf - t0)/dt);
%
% Find the indices of the measurement times nearest the initial and final
% times.
%
  i0 = find(tw >= (nd0-1)*24*3600 + t0,1,'first');
  iF = i0 + nt;
%
% Extract the time span.
%
  tspan = tw(i0:iF);
%
end

% Copyright 2016 Kevin J. Kircher. See bldg/license.txt for full license
% information.