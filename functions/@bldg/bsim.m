function [X,Q] = bsim(b,t,W,x0,varargin)
%BSIM simulates the simplest building.

% Syntax:
%	[X,Q] = bsim(b,t,W,x0) controls perfectly at setpoint Ts = x0(end).
%	[X,Q] = bsim(b,t,W,x0,options) uses specified ODE solver options.
%	X = bsim(b,t,W,x0,U) with size(U,2) == 1 uses specified control heat
%   flows Qc = U.
%	X = bsim(b,t,W,x0,U) with size(U,2) == 2 uses specified control air 
%   mass flow rates mdotc = U(:,1) at supply temperatures Tc = U(:,2).
%	X = bsim(b,t,W,x0,U,options) uses specified ODE solver options.
%
% Inputs:
%	b, a bldg object containing geometric, material, and usage parameters.
%	t, a nondecreasing array of nt simulation times in solar seconds.
%	W, a nt x 5-7 matrix containing the following disturbance histories:
%		W(:,1) = Tinf, the outdoor air temperature.
%		W(:,2) = Ih, the total horizontal irradiance.
%		W(:,3) = Ibn, the beam normal irradiance.
%		W(:,4) = Qp, the internal heat gained from people.
%		W(:,5) = Ql, the internal heat gained from lights.
%		W(:,6) = Qe, the internal heat gained from other equipment.
%	x0, an array of nx initial states.
%	U, the controls in either of two forms:
%		1) an nt x 1 vector of internal control heat gains Qc, or
%		2) an nt x 2 matrix containing:
%			U(:,1) = mdotc, the control air mass flow rate.
%			U(:,2) = Tc, the control air temperature.
%	options, a struct of ODE solver options defined by odeset.
%
% Outputs:
%	X, an nt x nx matrix containing the following state histories:
%		X(:,1) = T1, the indoor wall surface temeprature.
%		X(:,2:N-1) = [T2;...;TN-1], the internal wall node temperatures.
%		X(:,N) = TN, the outdoor wall surface temperature.
%		X(:,N+1) = Tg, the window glass temperature.
%		X(:,N+2) = Ta, the indoor air temperature.
%	Q, an array of nt internal heat gains from control systems, as
%		required to perfectly regulate air temperature at setpoint Ts.
%
% In the special case where bsim is called with the syntax
%   X = bsim(b,t,W,x0,U)
% with length(t) = 2, the inputs W and U may be vectors with length(W) = 6 
% and length(U) = 1 if U = Qc or length(U) = 2 if U = (mdotc, Tc). In this
% case, bsim acts as a discrete-time dynamics function. The output X is an
% column vector of the states at t(end).
%
% The input signals W and U are interpolated using a zero-order hold. For 
% all times tau satisfying t(i) < tau <= t(i+1), the control is U(i,:) and 
% the disturbance is W(i,:). When an ODE solver requires data before t(1)
% (or after t(end)), the data from t(1) (or t(end)) are used. 

%
% Check that the time vector is nondecreasing.
%
  if any(diff(t) <= 0)
    error('t must be increasing.')
  end
%
% Check that the syntax is correct. Determine whether perfect control is 
% required and whether ODE solver options are specified.
%
  nVariable = length(varargin);
  if nargout == 2                       % perfect control case
    perfectControl = 1;
    if nVariable == 0                   % no optios specified
      optionsSpecified = 0;
    elseif nVariable == 1               % options specified
      optionsSpecified = 1;
      options = varargin{1};
      if ~isstruct(options)
        error('Improper syntax: ODE options must be a struct defined by odeset.m.')
      end
    else
      error('Improper syntax: when called with 2 outputs, bsim requires 4 or 5 inputs.')
    end
  elseif nargout == 1                   % imperfect control case
    perfectControl = 0;
    if nVariable == 1                   % no options specified
      optionsSpecified = 0;
    elseif nVariable == 2               % options specified
      optionsSpecified = 1;
      options = varargin{2};
      if ~isstruct(options)
        error('Improper syntax: ODE options must be a struct defined by odeset.m.')
      end
    else
      error('Improper syntax: when called with 1 output, bsim requires 5 or 6 inputs.')
    end
    U = varargin{1};
    if isstruct(U)
      error('Improper syntax: control must be a matrix or array.')
    end
  else
    error('Improper syntax: bsim must be called with 1 or 2 outputs.')
  end	
%
% In the special case of the discrete-time dynamics call 
%   xf = bsim(b,t,w,x0,u), 
% check the dimensions of w and u and adjust them if necessary.
%
  nt = length(t); discreteTime = 0;
  if nt == 2 && ~perfectControl
    if isvector(W) && isvector(U)
      discreteTime = 1;
      if all(size(W) == [1,6])
        W = repmat(W,2,1);
      elseif all(size(W) == [6,1])
        W = repmat(W',2,1);
      elseif any(size(W) ~= [2,6])
        error('For length(t) = 2, the disturbance W must be either 2 x 6 or a 6-vector.')
      end
      if all(size(U) == [1,1]) || all(size(U) == [1,2])
        U = repmat(U,2,1);
      elseif all(size(U) == [1,2])
        U = repmat(U',2,1);
      elseif ~(all(size(U) == [2,1]) || all(size(U) == [2,2]))
        error('For length(t) = 2, the control U must be 1 x 1, 1 x 2, or 2 x 1.')
      end
    end
  end
%
% Make sure the disturbance and control dimensions are correct.
%
  if size(W,1) ~= nt
    error('For length(t) = nt > 2, the disturbance W must be nt x 6.')
  end
  if ~perfectControl && (size(U,1) ~= nt || ~ismember(size(U,2),[1,2]))
    error('For length(t) = nt > 2, the control U must be nt x 1 or nt x 2.')
  end
%
% Extract relevant geometry and material properties.
%
  l = b.l;                              % wall thickness
  A = b.A;                              % glass/wall surface area
  a = b.a;                              % wall diffusivity
  k = b.k;                              % wall conductivity
  Cg = b.Cg;                            % glass capacitance
  Ca = b.Ca;                            % air capacitance
  ca = b.ca;                            % air specific heat at constant pressure
%
% Define numerical parameters and finite difference matrix.
%
  N = length(x0) - 2;
  dx = l/(N-1);                         % mesh width
  r = a/dx^2;                           % mesh ratio
  mainDiag = -2*r*ones(N,1);
  superDiag = zeros(N,1); superDiag(2) = 2*r; superDiag(3:end) = r;
  subDiag = zeros(N,1); subDiag(end-1) = 2*r; subDiag(1:end-2) = r;
  Afd = spdiags([subDiag,mainDiag,superDiag],-1:1,N,N);
%
% Define sparsity pattern of Jacobian matrix.
%
  Df = zeros(N+2,N+2);
  Df(1:N,1:N) = full(Afd); 
  Df(Df ~= 0) = 1;                      % wall nodes to themselves
  Df(1,N+1) = 1;                        % inner wall surface to glass
  Df(1,N+2) = 1;                        % inner wall surface to air
  Df(N+1,1) = 1;                        % glass to inner wall surface
  Df(N+1,N+1) = 1;                      % glass to glass
  Df(N+1,N+2) = 1;                      % glass to air
  if ~perfectControl
    Df(N+2,1) = 1;                      % air to inner wall surface
    Df(N+2,N+1) = 1;                    % air to glass
    Df(N+2,N+2) = 1;                    % air to air
  end
%
% Add the Jacobian sparsity pattern to the ODE solver options and flag the 
% dynamics function f(t,x) as vectorized.
%
  if optionsSpecified
    options = odeset(options,'JPattern',Df, 'Vectorized','on');
  else
    options = odeset('JPattern',Df, 'Vectorized','on');
  end
%
% Numerically integrate the building state.
%
  if perfectControl
    [~,X] = ode15s(@f,t,x0,options);		% ode15s can handle DAEs
  else
    [~,X] = ode23tb(@f,t,x0,options);		% ode23tb cannot
  end
  if nt == 2
    if discreteTime
      X = X(end,:)';
    else
      X = [X(1,:);X(end,:)]; 
    end
  end
%
% If necessary, compute controls required for perfect regulation.
%
  if perfectControl
    Q = getLoads(b,X,W,t);
  end	
%
%
%----------------------- Nested dynamics function. ------------------------
%
% 
  function xdot = f(tNum,x)
  %F defines the dynamics of the simplest building.

  % Inputs: 
  %		tNum, the numerical solution time in seconds.
  %		x, the building state, contains
  %			x(1), the temperature of the inner wall surface.
  %			x(2), ..., x(N-1), the temperatures inside the wall.
  %			x(N), the temperature of the outer wall surface.
  %			x(N+1), the glass temperature.
  %			x(N+2), the indoor air temperature.
  %
  % Output:
  %		xdot, the time derivative of the state.

  %
  % Get the weather and internal gains.
  %
    [Tinf,Id,Ibp,Ibg,thg] = getWeather(b,W,t,tNum);
    if perfectControl
      [Qp,Ql,Qe] = getGains(W,t,tNum);
    else
      [Qp,Ql,Qe,Qc] = getGains(W,t,tNum,ca,x,U);
    end	  
  %
  % Compute the surface fluxes, including infiltration.
  %
    if perfectControl
      [qgm,qgp,qm,qp] = getFluxes(b,x,Tinf,Id,Ibp,Ibg,thg,Qp,Ql,Qe);
    else
      [qgm,qgp,qm,qp,qam,qap,Qa] = getFluxes(b,x,Tinf,Id,Ibp,Ibg,thg,Qp,Ql,Qe,Qc);
    end
  %
  % Compute the wall temperature derivatives.
  %
    [nx,ncol] = size(x);
    xdot = zeros(nx,ncol);
    xdot(1:N,:) = Afd*x(1:N,:);
    xdot(1,:) = xdot(1,:) + 2*a*qm'/(k*dx);
    xdot(N,:) = xdot(N,:) - 2*a*qp'/(k*dx);
  %
  % Compute the glass and air temperature derivatives.
  %
    xdot(N+1,:) = A*(qgm-qgp)'/Cg;				% glass
    if perfectControl                     % air
      xdot(N+2,:) = 0;
    else
      xdot(N+2,:) = (A*(qam-qap)+Qa)'/Ca;
    end
  %
  end
%
end

% Copyright 2016 Kevin J. Kircher. See bldg/license.txt for full license
% information.