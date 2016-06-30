function [b,weather] = importWeather(b,fileName,dt)
%IMPORTWEATHER imports and interpolates weather data from a TMY3 file.

% Inputs:
%	b, a bldg object containing building parameters.
%	fileName, a string containing the the name of the weather file.
%	dt, the sample time.
%
% Outputs:
%	b, the input bldg object with the following property modified:
%		phi, the latitude in radians.
%	weather, a struct containing:
%		tw, the measurement solar times in seconds.
%		Tinf, the outdoor air dry bulb temperatures.
%		Ih, the total horizontal irradiance.
%		Ibn, the beam normal irradiance.

%
% Check that the weather file is on the path.
%
  if ~exist(fileName,'file')
    error('The weather file must be on your path or in the working directory.')
  end
%
% Make sure the sample time is no more than one hour.
%
  if dt > 60*3600
    error('Sample times over one hour are not currently supported.')
  end
%
% Read the latitude and longitude from the first line.
%
  fid = fopen(fileName);
  angles = textscan(fid,'%*s %*s %*s %*s %f %f %*s \n','Delimiter',',');
  b.phi = pi*angles{1}/180;                  % latitude
  psi = pi*abs(angles{2})/180;                  % longitude
  fclose(fid);
%
% Read the solar and temperature data.
%
  fid = fopen(fileName);
  format = strcat('%*s %*s %*f %*f %f %*f %*f %f',repmat(' %*s',1,23),' %f',...
	repmat(' %*s',1,39));
  C = textscan(fid,format,'Delimiter',',','HeaderLines',2);
  fclose(fid);
  Ih = C{1};                                    % total horizontal irradiance
  Ibn = C{2};                                   % beam normal irradiance
  Tinf = 273 + C{3};                            % dry bulb temperature
%
% Determine the solar time of the measurements, in seconds.
%
  psiDeg = 180*psi/pi;
  psiStd = round(psiDeg/15)*15;
  dPsi = psiStd - psiDeg;
  nh = length(Ih);
  hc = (0:nh-1)';                               % clock hours
  tc = hc*3600;                                 % clock seconds
  nd = floor(tc/(24*3600)) + 1;	nd(end) = 365;	% day numbers
  lam = 2*pi*(nd - 1)/365;                      % day angles
  tw = tc + ( 4*dPsi + 229.2*(0.000075 + 0.001868*cos(lam) ...
    - 0.032077*sin(lam) - 0.014615*cos(2*lam) - 0.04089*sin(2*lam)))/60;
  tw = tw - tw(1);
%
% Interpolate the data and store them in the output struct.
%
  Imin = 7;                                     % irradiance threshold
  method = 'spline';                            % interpolation method
  if dt < 60*3600
    nt = floor(nh*3600/dt);
    temp_tw = (tw(1):dt:tw(end))';
    deficit = nt - length(temp_tw);
    if deficit > 0
      weather.tw = [temp_tw;temp_tw(end)+dt*(1:deficit)'];
    elseif deficit < 0
      weather.tw = temp_tw(1:nt);
    else
      weather.tw = temp_tw;
    end
    weather.tw = weather.tw(1:nt);
    weather.Tinf = interp1(tw,Tinf,weather.tw,method);
    weather.Ih = max(0,interp1(tw,Ih,weather.tw,method));
    weather.Ih(weather.Ih < Imin) = 0;					% interpolation clean-up
    weather.Ibn = max(0,interp1(tw,Ibn,weather.tw,method));
    weather.Ibn(weather.Ibn < Imin) = 0;				% interpolation clean-up
  else
    weather.tw = tw;
    weather.Tinf = Tinf;
    weather.Ih = Ih;
    weather.Ibn = Ibn;
  end
%
end

% Copyright 2016 Kevin J. Kircher. See bldg/license.txt for full license
% information.