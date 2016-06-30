function gains = generateGains(b,t,varargin)
%GENERATEGAINS generates stochastic internal gains for the simplest building.

% Syntax:
% gains = generateGains(b,t) generates gains for a residential building.
% gains = generateGains(b,t,buildingType) generates gains for either a
%   residential or commercial building, depending on buildingType.
%
% Inputs:
%	b, a bldg object containing:
%		eta, the lighting efficiency.
%	t, the solar time span in seconds.
% buildingType, a string with either the value 'residential' or
%   'commercial'.
%
% Output:
%	gains, a struct containing:
%		Qp, a column of body heat flows.
%		Ql, a column of lighting heat flows.
%		Qe, a column of equipment heat flows.

%
% Parse the inputs to decide whether to generate gains for a residential or
% commercial building.
%
  nVariable = length(varargin);
  if nVariable == 0
    bldgType = 'residential';       % default: residential
  elseif nVariable == 1
    thirdArgument = varargin{1};
    if ~ischar(thirdArgument)
      error('buildingType must be a character array.')
    else
      if strcmpi(thirdArgument,'residential')
        bldgType = 'residential';
      elseif strcmpi(thirdArgument,'commercial')
        bldgType = 'commercial';
      else
        error('buildingType must be either residential or commercial.')
      end
    end
  else
    error('generateGains must be called with either 2 or 3 arguments.')
  end
    
%
% Extract relevant fields from input struct.
%
  eta = b.eta;
%

%% Tunable parameters
%
% People.
%
  bp = 120;                         % body heat per person
  footprint = b.la^2;               % building footpring (m^2)
  nFloors = floor(b.A/b.la/3);      % number of floors
  areaPerPerson = 100;
  np = max(1,round(footprint*nFloors...
    /areaPerPerson));               % maximum number of occupants
%
% Lights. (Default: CFLs, typical office building capacity.)
%
  lightPowerIntensity = 10;          % W/m^2
  Pl0 = lightPowerIntensity*...
    footprint*nFloors;              % peak lighting electric power, W
  night = [0,7];                    % quiet hours
%
% Equipment.
%
  equipPowerIntensity = 10;         % W/m^2
  Qe0 = equipPowerIntensity*...
    footprint*nFloors;              % constant (offset), W
  be = 300;                         % equipment heat per occupant (slope)
% 
% Autocorrelation
%
  rp = 0.9;                         % 15min: 0.9; 5min:
  re = -0.4;                        % 15min: 0.4; 5min:
%

%% Gain generation
%
% Initialize.
%
  nt = length(t);
  p = zeros(nt,1);
%
% Find the set of weekdays. Take January 1 to be Monday.
%
  nd = ceil(t/(24*3600));           % day number
  nd(1) = 1; nd(end) = 365;
  day = mod(nd-1,7)+1;              % day of week
  wkd = zeros(nt,1);                % weekday = 0
  wkd(day > 5) = 1;                 % weekend = 1
%
% Loop through the days, defining occupancy.
%
  td = linspace(0,24*3600,nt/365);
  for id = 1 : 365
    occ = zeros(nt/365,1);
	  it = (id-1)*nt/365 + 1;
    if wkd(it)                        % weekend pattern
      dep = (9+randi([0,9],1))*3600;
      arr = (19+randi([0,5],1))*3600;
    else                              % weekday pattern
      dep = (7+randi([0,2],1))*3600;
      arr = (16+randi([0,4],1))*3600;
    end
    occ(td < dep | td > arr) = 1;
    p(it:it+nt/365-1) = occ;
  end
%
% If the building type is commercial, flip occupancy for a workplace.
%
  if strcmpi(bldgType,'commercial')
    pnew = 0*p;
    pnew(p == 0) = 1;
    p = pnew;
  end
%
% Initialize the light and equipment heat flows.
%
  Ql = zeros(nt,1);
  Qe = zeros(nt,1);
%
% Compute the gains.
%
  small = 1e-8;
  Qp = bp*np*p;
  h = mod(t,24*3600)/3600;					% hours since midnight
  Ql(abs(p) > small & (h < night(1) | h > night(end))) = (1-eta)*Pl0;
  Qe(h < night(1) | h > night(end)) = be*np*p(h < night(1) | h > night(end));
  Qe = Qe + Qe0;
%
% Corrupt the heat flows with autocorrelated noise. Store the data in the
% output struct.
%
  sigp = 0.10*bp;
  sigl = 0.05*mean(Ql);
  sige = 0.10*Qe0;
  tempp = Qp + sigp*filter(1,[1,-rp],randn(nt,1)); tempp(p == 0) = 0;
  gains.Qp = max(0,tempp);
  templ = Ql + sigl*randn(nt,1); templ(Ql == 0) = 0;
  gains.Ql = max(0,templ);
  gains.Qe = max(0,Qe + sige*filter(1,[1,-re],randn(nt,1)));
%
end

% Copyright 2016 Kevin J. Kircher. See bldg/license.txt for full license
% information.