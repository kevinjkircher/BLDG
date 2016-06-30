function Qc = getLoads(b,X,W,t)
%GETLOADS gets the heating/cooling load of the simplest building.

% Inputs:
%	b, a bldg object.
%	X, an nt x nx matrix containing the following state histories:
%		X(:,1) = T1, the indoor wall surface temeprature.
%		X(:,2:N-1) = [T2;...;TN-1], the internal wall node temperatures.
%		X(:,N) = TN, the outdoor wall surface temperature.
%		X(:,N+1) = Tg, the window glass temperature.
%		X(:,N+2) = Ta, the indoor air temperature.
%	W, a nt x 5 matrix containing the following disturbance histories:
%		W(:,1) = Tinf, the outdoor air temperature.
%		W(:,2) = Ih, the total horizontal irradiance.
%		W(:,3) = Ih, the beam normal irradiance.
%		W(:,4) = Qp, the internal heat gained from people.
%		W(:,5) = Ql, the internal heat gained from lights.
%		W(:,6) = Qe, the internal heat gained from other equipment.
%	t, a monotonic array of nt simulation times in solar seconds.
%
% Output:
%	Qc, the heat flows from control systems.

%
% Extract parameters, weather, and internal gains.
%
  A = b.A;                            % glass/wall surface area
  ca = b.ca;                          % air specific heat
  mdot = b.mdot;                      % infiltration mass flow rate
  zc = b.zc;                          % control convective fraction
  zp = b.zp;                          % people convective fraction
  zl = b.zl;                          % lighting convective fraction
  ze = b.ze;                          % equipment convective fraction
  [Tinf,~,~,~,~] = getWeather(b,W,t,t);  
  [Qp,Ql,Qe] = getGains(W,t,t);
%
% Extract temperatures.
%
  Ta = X(:,end);
  Tm = X(:,1);                        % inner wall surface
  Tg = X(:,end-1);                    % inner glass surface
%
% Compute film coefficients and surface fluxes.
%
  [hm,hgp] = getFilmCoefficients(Ta,Tg,Tm);
  qam = hgp.*(Tg - Ta);
  qap = hm.*(Ta - Tm);
%
% Compute loads.
%
  Qc = (A*(qap-qam) - (mdot.*ca.*(Tinf-Ta) + zp*Qp + zl*Ql + ze*Qe))/zc;
%
end

% Copyright 2016 Kevin J. Kircher. See bldg/license.txt for full license
% information.