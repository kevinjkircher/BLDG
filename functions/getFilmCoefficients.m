function [hm,hgp,hgm,hp] = getFilmCoefficients(Ta,Tg,Tm,Tp,Tinf)
%GETFILMCOEFFICIENTS defines the film (convection) coefficients for the 
%simplest building.

% Syntax:
%	[hm,hgp] = getFilmCoefficients(Ta,Tg,Tm)
%	[hm,hgp,hgm,hp] = getFilmCoefficients(Ta,Tg,Tm,Tp,Tinf)
%
% Inputs:
%	Ta, the indoor air temperature.
%	Tg, the glass temperature.
%	Tm, the inner wall surface temperature.
%	Tp, the outer wall surface temperature.
%	Tinf, the outdoor air temperature.
%
% Outputs:
%	hgm, the film coefficient at the left (outdoor) glass surface.
%	hgp, the film coefficient at the right glass surface.
%	hm, the film coefficient at the left (indoor) wall surface.
%	hp, the film coefficient at the right wall surface.

%
% Define a general film coefficient model.
%
  h = @(T1,T2,C1,C2) C1*abs(T1-T2).^C2;
%
% Set the model parameters for each surface using an appropriate empirical
% model from the literature.
%
  C1_ashrae = 1.31; C2_ashrae = 1/3;            % ASHRAE model
  C1_khalifa = 1.98; C2_khalifa = 0.32;         % Khalifa 1989 model
%   C1_trnsys = 1.5; C2_trnsys = 0.25;            % TRNSYS model
%
% Compute the film coefficients.
%
  hm = h(Tm,Ta,C1_khalifa,C2_khalifa);
  hgp = h(Tg,Ta,C1_khalifa,C2_khalifa);
  if nargin > 3
    hgm = h(Tg,Tinf,C1_ashrae,C2_ashrae);
    hp = h(Tp,Tinf,C1_ashrae,C2_ashrae);
  end
%
end

% Copyright 2016 Kevin J. Kircher. See bldg/license.txt for full license
% information.