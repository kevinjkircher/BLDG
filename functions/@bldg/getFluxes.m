function [qgm,qgp,qm,qp,qam,qap,Qa] = getFluxes(b,x,Tinf,Id,Ibp,Ibg,thg,Qp,Ql,Qe,Qc)
%GETFLUXES gets the surface fluxes for the simplest building.

% Syntax:
%	[qgm,qgp,qm,qp] = getFluxes(b,x,Tinf,Id,Ibp,Ibg,thg,Qp,Ql,Qe)
%		computes fluxes assuming perfect control.
%	[qgm,qgp,qm,qp,qam,qap,Qa] = getFluxes(b,x,Tinf,Id,Ibp,Ibg,thg,Qp,Ql,Qe,Qc)
%		computes fluxes under user-specified control.
%
% Inputs:
%	b, a bldg object.
%	x, the building state.
%	Tinf, the outdoor air temperature.
%	Id, the diffuse solar irradiance, assumed isotropic.
%	Ibp, the beam solar irradiance on the outdoor wall surface.
%	Ibg, the beam solar irradiance on the window.
%	thg, the beam angle of incidence on the window.
%	Qp, the heat flow from people.
%	Ql, the heat flow from lights.
%	Qe, the heat flow from equipment.
%
% Outputs:
%	qgm, the left (outdoor) glass surface flux.
%	qgp, the right glass surface flux.
%	qam, the left air surface flux.
%	qap, the right air surface flux.
%	qm, the left wall surface flux.
%	qp, the right (outdoor) wall surface flux.
%	Qe, the total heat flow into the air from internal gains.

%
% Infer whether perfect or user-specified control is in use.
%
  perfectControl = 1;
  if nargin > 10, perfectControl = 0; end
%
% Extract parameters.
%
  N = size(x,1) - 2;
  ca = b.ca;                          % air specific heat
  sig = b.sig;                        % Stefan-Boltzmann constant
  lg = b.lg;                          % glass thickness
  A = b.A;                            % glass/wall surface area
  eps = b.eps;                        % wall longwave emissivity
  as = b.as;                          % wall shortwave absorptivity
  epsg = b.epsg;                      % glass longwave emissivity
  ng = b.ng;                          % glass index of refraction
  mug = b.mug;                        % glass attenuation coefficient
  agBar = b.agBar;                    % glass diffuse absorptivity
  taugBar = b.taugBar;                % glass diffuse transmissivity
  mdot = b.mdot;                      % infiltration mass flow rate
  zc = b.zc;                          % control convective fraction
  zp = b.zp;                          % people convective fraction
  zl = b.zl;                          % lighting convective fraction
  ze = b.ze;                          % equipment convective fraction
  eta = b.eta;                        % lighting efficiency
%
% Extract states. Store as columns.
%
  Tm = x(1,:)';                       % inner wall surface
  Tp = x(N,:)';                       % outer wall surface
  Tg = x(N+1,:)';                     % glass
  Ta = x(N+2,:)';                     % indoor air
%
% Compute film coefficients and glass beam absorptivity/emissivity.
%
  [hm,hgp,hgm,hp] = getFilmCoefficients(Ta,Tg,Tm,Tp,Tinf);
  [agthg,taugthg] = getAbsTrans(thg,lg,ng,mug);
%
% Compute fluxes from window and wall into air. Solve for the control heat
% in the perfect control case.
%
  qam = hgp.*(Tg - Ta);
  qap = hm.*(Ta - Tm);
  if perfectControl
    Qc = (A*(qap-qam) - (mdot.*ca.*(Tinf-Ta) + zp*Qp + zl*Ql + ze*Qe))/zc;
  end
%
% Compute external surface fluxes.
%
  qp = hp.*(Tp - Tinf) + sig*eps*(Tp.^4 - Tinf.^4) - as.*(Id + Ibp);
  qgm = hgm.*(Tinf - Tg) + sig*epsg*(Tinf.^4 - Tg.^4) + agthg.*Ibg + agBar.*Id;
%
% Compute internal shortwave fluxes.
%
  L = eta*Ql/(A*(1-eta));
  rhogBar = 1 - agBar - taugBar;			% glass diffuse reflectivity
  rhos = 1 - as;                      % wall shortwave reflectivity
  qmShort = as*(taugthg.*Ibg + taugBar*Id + (1 + rhogBar)*L)...
    / (1 - rhos*rhogBar);
  qgpShort = - agBar*(rhos*(taugthg.*Ibg + taugBar*Id) + (1 + rhos)*L)...
    / (1 - rhos*rhogBar);
%
% Compute internal longwave fluxes.
%
  sIn = (0.5/A)*((1-zc)*Qc + (1-zp)*Qp + (1-zl)*Ql + (1-ze)*Qe);
  qmIn = (2-epsg)*sIn / (1+epsg*(1/eps-1));
  qgpIn = - (2-eps)*sIn / (1+eps*(1/epsg-1));
  qgwNet = sig*(Tg.^4 - Tm.^4) / (1/epsg + 1/eps - 1);
  qmLong = qmIn + qgwNet;
  qgpLong = qgpIn + qgwNet;
%
% Compute total internal surface fluxes.
%
  qgp = qam + qgpShort + qgpLong;
  qm = qap + qmShort + qmLong;
%
% Compute total heat flow into the air from internal gains.
%
  if ~perfectControl
    Qa = mdot.*ca.*(Tinf - Ta) + zc*Qc + zp*Qp + zl*Ql + ze*Qe;
  end
%
end

% Copyright 2016 Kevin J. Kircher. See bldg/license.txt for full license
% information.