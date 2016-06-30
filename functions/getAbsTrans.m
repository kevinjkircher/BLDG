function [agthg,taugthg] = getAbsTrans(thg,lg,ng,mug)
%GETABSTRANS calculates the absorptivity and transmissivity of glass with
%respect to beam solar irradiance.

% Inputs: 
%	tgh, the beam angle of incidence.
%	lg, the glass thickness.
%	ng, the glass index of refraction.
%	mug, the glass attenuation coefficient.
%
% Output:
%	agthg, the beam absorptivity.
%	taugthg, the beam transmissivity.

%
% If any angles of incidence are zero, add a little bit to them.
%
  small = 1e-8;
  thg(thg < small) = small;
%
% Compute the reflected and absorbed fractions.
%
  thr = asin(sin(thg)/ng);					% angle of refraction
  rthg = 0.5*(tan(thg - thr).^2 ./ tan(thg + thr).^2 ...
	+ sin(thg - thr).^2 ./ sin(thg + thr).^2);% reflected fraction
  athg = 1 - exp(-mug*lg ./ cos(thr));		% absorbed fraction
%
% Compute the absorptivity and transmissivity.
%
  agthg = athg.*(1 - rthg) ./ (1 - rthg.*(1 - athg));
  taugthg = (1 - rthg).^2.*(1 - athg) ./ (1 - rthg.^2.*(1 - athg).^2);
%
end

% Copyright 2016 Kevin J. Kircher. See bldg/license.txt for full license
% information.