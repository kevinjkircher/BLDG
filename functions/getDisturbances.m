function W = getDisturbances(weather,gains,t)
%GETDISTURBANCES builds a matrix of distrubances from the weather, gains,
%and simulation time span, for use with bsim.

% Inputs:
%	weather, a struct containing:
%		tw, the measurement times in solar seconds.
%		Tinf, a column of dry bulb air temperatures.
%		Ih, a column of total horizontal irradiances.
%		Ibn, a column of beam normal irradiances.
%	gains, a struct containing:
%		Qp, a column of body heat flows.
%		Ql, a column of lighting heat flows.
%		Qe, a column of equipment heat flows.
%	t, the simulation solar time span.
%
% Output:
%	W, a matrix of disturbances. Each column contains the values of a
%	disturbance at each time in t.
%		W(1,:) = Tinf
%		W(2,:) = Ih
%		W(3,:) = Ibn
%		W(4,:) = Qp
%		W(5,:) = Ql
%		W(6,:) = Qe

%
% Make sure t is a column vector.
%
  if size(t,1) == 1, t = t'; end
%
% Build the disturbance matrix.
%
  [~,indices] = histc(t,weather.tw);		% indices of nearest measurements
  nd = 5;									% number of disturbances
  nt = length(t);							% number of measurement times
  W = zeros(nt,nd);
  W(:,1) = weather.Tinf(indices);
  W(:,2) = weather.Ih(indices);
  W(:,3) = weather.Ibn(indices);
  W(:,4) = gains.Qp(indices);
  W(:,5) = gains.Ql(indices);
  W(:,6) = gains.Qe(indices);
%
end

% Copyright 2016 Kevin J. Kircher. See bldg/license.txt for full license
% information.