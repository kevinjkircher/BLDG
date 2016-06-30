function [Qp,Ql,Qe,Qc] = getGains(W,t,tNum,ca,x,U)
%GETGAINS gets the internal gains for the simplest building.

% Syntax:
%	[Qp,Ql,Qe] = getGains(W,t,tNum)
%	[Qp,Ql,Qe,Qc] = getGains(W,t,tNum,ca,x,U)
%
% Inputs:
%	W, a nt x 5 matrix containing the following disturbance histories:
%		W(:,1) = Tinf, the outdoor air temperature.
%		W(:,2) = Ih, the total horizontal irradiance.
%		W(:,3) = Ibn, the beam normal irradiance.
%		W(:,4) = Qp, the internal heat gained from people.
%		W(:,5) = Ql, the internal heat gained from lights.
%		W(:,6) = Qe, the internal heat gained from other equipment.
%	t, the full array of the simulation times.
%	tNum, an array of times of interest to the ODE solver.
%	ca, the specific heat of air.
%	x, the building state.
%	U, the controls in either of two forms:
%		1) an nt x 1 vector of internal control heat gains Qc, or
%		2) an nt x 2 matrix containing:
%			U(:,1) = mdotc, the control air mass flow rate.
%			U(:,2) = Tc, the control air temperature.
%
% Outputs:
%	Qp, the heat flow from people.
%	Ql, the heat flow from lights.
%	Qe, the heat flow from equipment.
%	Qc, the heat flow from control systems.

%
% Infer whether perfect or user-specified control is in use.
%
  perfectControl = 0;
  if nargin < 4, perfectControl = 1; end
%
% Make sure input vectors are columns.
%
  if size(t,1) == 1, t = t'; end
  if size(tNum,1) == 1, tNum = tNum'; end
%
% For any times before (after) the first (last) simulation time, use the 
% first (last) time.
%
  tNum(tNum < t(1)) = t(1);
  tNum(tNum > t(end)) = t(end);
%
% Get the nearest data on gains from people, lights and equipment.
%
  [~,indices] = histc(tNum,t);            % indices of nearest data
  Qp = W(indices,4);
  Ql = W(indices,5);
  Qe = W(indices,6);
%
% Get the controls.
%
  if ~perfectControl
    if size(U,2) == 1
      Qc = U(indices);
    elseif size(U,2) == 2
      N = size(x,1) - 2;                  % number of nodes
      Ta = x(N+2,:)';                     % room air temperature
      mdotc = U(indices,1);               % control mass flow rate
      Tc = U(indices,2);                  % control air temperature
      Qc = mdotc.*ca.*(Tc - Ta);
    else
      error('Size(U,2) must be either 1 (if U = Qc) or 2 (if U = [mdotc,Tc]).')
    end
  end
%
end

% Copyright 2016 Kevin J. Kircher. See bldg/license.txt for full license
% information.