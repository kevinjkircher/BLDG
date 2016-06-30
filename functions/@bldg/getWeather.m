function [Tinf,Id,Ibp,Ibg,thg] = getWeather(b,W,t,tNum)
%GETWEATHER gets weather for the simplest building.

% Inputs: 
%	b, a bldg object containing the following fields:
%		gam, the wall surface azimuth angle in radians (0=south, pi/2=west).
%		phi, the latitude in radians.
%	W, a nt x 5 matrix containing the following disturbance histories:
%		W(:,1) = Tinf, the outdoor air temperature.
%		W(:,2) = Ih, the total horizontal irradiance.
%		W(:,3) = Ih, the beam normal irradiance.
%	t, the full array of the simulation times.
%	tNum, an array of times of interest to the ODE solver.
%
% Outputs:
%	Tinf, the outdoor air temperature.
%	Id, the diffuse solar irradiance.
%	Ibp, the beam solar irradiance on the wall.
%	Ibg, the beam solar irradiance on the window.
%	thg, the beam angle of incidence on the window.

%
% Make sure time vectors are columns.
%
  if size(t,1) == 1, t = t'; end
  if size(tNum,1) == 1, tNum = tNum'; end
%
% Get the wall azimuth and latitude. Compute the day numbers and angles.
%
  gam = b.gam;                                        % wall azimuth
  phi = b.phi;                                        % latitude
  nd = floor(tNum/(24*3600)) + 1;                     % day numbers
  lam = 2*pi*(nd - 1)/365;                            % day angles
%
% For any times before (after) the earliest (latest) simulation time, use 
% the earliest (latest) time.
%
  tNum(tNum < t(1)) = t(1);
  tNum(tNum > t(end)) = t(end);
%
% Get the nearest temperature and irradiance data.
%
  [~,indices] = histc(tNum,t);                        % indices of nearest data
  Tinf = W(indices,1);
  Ih = W(indices,2);
  Ibn = W(indices,3);
%
% Compute the solar zenith angle.
%
  del = 0.409*sin(lam + 4.9);                         % declination
  s = mod(tNum,24*3600);                              % seconds since midnight
  om = 0.262*(s/3600 - 12);                           % hour angle
  cosxi = cos(phi).*cos(del).*cos(om) ...             % zenith angle
    + sin(phi).*sin(del);					
%
% Split out the horizontal beam and diffuse components.
%
  small = cos(87*pi/180);                             % threshold avoids div_zero
  Ibh = Ibn.*cosxi;                                   % beam
  Ibh(cosxi < small) = 0;
  Idh = Ih - Ibh;                                     % diffuse
%
% Compute the beam and diffuse components on the wall and window.
%
  Id = Idh;                                           % assume isotropic
  costh = sin(phi).*cos(gam).*cos(del).*cos(om) ...		% wall angle of incidence
    - cos(phi).*cos(gam).*sin(del) + sin(gam).*cos(del).*sin(om);
  costhg = -costh;                                    % window angle of incidence
  thg = acos(costhg);                                 % ""
  Ibp = Ibh.*costh./cosxi;
  Ibg = Ibh.*costhg./cosxi;
  Ibp(cosxi < small | costh < small) = 0;             % dark or occluded times
  Ibg(cosxi < small | costhg < small) = 0;            % dark or occluded times
%
% Make sure the irradiances are positive.
%
  Id(Id < 0) = 0;
  Ibp(Ibp < 0) = 0;
  Ibg(Ibg < 0) = 0;
%
end

% Copyright 2016 Kevin J. Kircher. See bldg/license.txt for full license
% information.