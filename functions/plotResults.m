function plotResults(t,X,W,U,Tmin,Tmax)
%PLOTRESULTS plots the states and controls of the simplest building.

% Syntax:
%	plotResults(t,X,W,U) plots temperatures and controls vs. time.
%	plotResults(t,X,W,U,Tmin,Tmax) adds lines at T=Tmin and T=Tmax.
%
% Inputs:
%	t, a nondecreasing array of nt simulation times in solar seconds.
%	X, an nt x nx matrix containing the following state histories:
%		X(:,1) = T1, the indoor wall surface temeprature.
%		X(:,2:N-1) = [T2;...;TN-1], the internal wall node temperatures.
%		X(:,N) = TN, the outdoor wall surface temperature.
%		X(:,N+1) = Tg, the window glass temperature.
%		X(:,N+2) = Ta, the indoor air temperature.
%	W, a nt x 5 matrix containing the following disturbance histories:
%		W(:,1) = Tinf, the outdoor air temperature.
%		W(:,2) = Ih, the total horizontal irradiance.
%		W(:,3) = Ibn, the beam normal irradiance.
%		W(:,4) = Qp, the internal heat gained from people.
%		W(:,5) = Ql, the internal heat gained from lights.
%		W(:,6) = Qe, the internal heat gained from other equipment.
%	U, the controls in either of two forms:
%		1) an nt x 1 vector of internal control heat gains Qc, or
%		2) an nt x 2 matrix containing:
%			U(:,1) = mdotc, the control air mass flow rate.
%			U(:,2) = Tc, the control air temperature.
%	Tmin, the lower limit of the temperature deadband. in Kelvin.
%	Tmax, the upper limit of the temperature deadband, in Kelvin.

%
% Get the timing.
%
  N = size(X,2) - 2;
  nd0 = floor(t(1)/(24*3600)) + 1;					% initial day number
  s = t - (nd0-1)*24*3600;							% seconds since midnight
  h = s/3600;										% hours since midnight
  nh = length(h);
  hMax = ceil(h(end));
  if ceil(hMax) > 3*4*7*24, unit = 'Month'; x = h/24/7/4; xMax = hMax/24/30.4;
  elseif ceil(hMax) > 3*7*24, unit = 'Week'; x = h/24/7; xMax = hMax/24/7;
  elseif ceil(hMax) > 3*24, unit = 'Day'; x = h/24; xMax = hMax/24;
  else unit = 'Hour'; x = h; xMax = hMax; end
%
% Set the axis limits and ticks and the font size for titles and labels.
%
  xLimits = [round(x(1)),xMax];
  if strcmp(unit,'Month') || strcmp(unit,'Week') || strcmp(unit,'Day')
    xTicks = xLimits(1):floor(xLimits(end));
  elseif strcmp(unit,'Hour')
    xTicks = xLimits(1):6:ceil(xLimits(end));
  end
  fs = 16;											% label font size
%
% Get the default interpreter. Temporarily change it to LateX.
%
  original_interpreter = get(0,'DefaultTextInterpreter');
  set(0,'DefaultTextInterpreter','LaTeX');
%
% Plot the temperatures.
%
  clf, subplot(3,1,1:2), hold on
  XC = X - 273; TinfC = W(:,1) - 273;
  plot(x,XC(:,1),'r', x,XC(:,N),'m', x,XC(:,N+1),'g', x,XC(:,N+2),'k', ...
    x,TinfC,'b','linewidth',2)
  box on, xlim(xLimits), ylabel('Temperature (C)','FontSize',fs)
  yMin = floor(min(min([XC(:,[1 N:N+2]) TinfC]))/5)*5;
  yMax = ceil(max(max([XC(:,[1 N:N+2]) TinfC]))/5)*5;
  yLimits = [yMin yMax];
  ylim(yLimits)
  leg = legend('$T_-$','$T_+$','$T_g$','$T_a$','$T_\infty$',...
    'Location','NorthOutside','Orientation','Horizontal');
  set(leg,'FontSize',fs+2,'Interpreter','LaTeX')
  set(gca,'XTick',xTicks);
  set(gca,'XTickMode','Manual'), set(gca,'YTickMode','Manual')
%
% If required, add lines at T=Tmin and T=Tmax.
%
  if nargin > 4
    Tmins = (Tmin - 273)*ones(nh,1);
    Tmaxes = (Tmax - 273)*ones(nh,1);
    hold on, plot(x,Tmins,'k--', x,Tmaxes,'k--')
  end
%
% Add a subplot with the control(s).
%
  subplot(3,1,3)
  if size(U,2) == 1, forcedAir = 0;
  elseif size(U,2) == 2, forcedAir = 1;
  else error('size(u,2) must be either 1 or 2.')
  end
  if ~forcedAir											% u = Qc
    stairs(x,U/1000,'k','linewidth',2)
    xlim(xLimits)
    ylabel('$Q_c$ (kW)','FontSize',fs)
    set(gca,'XTick',xTicks);
    set(gca,'XTickMode','Manual'), set(gca,'YTickMode','Manual')
  else													% u = [mdotc, Tc]
    mdotc = U(:,1); Tc = U(:,2) - 273;
    [hAx,h1,h2] = plotyy(x,mdotc, x,Tc);
    set(h1,'linewidth',2); set(h2,'linewidth',2);
    ylabel(hAx(1),'$\dot{m}_c$ (kg/s)','FontSize',fs)	% left axis
    ylabel(hAx(2),'$T_c$ (C)','FontSize',fs)			% right axis
    xlim(hAx(1),xLimits), xlim(hAx(2),xLimits)
    set(hAx(1),'XTick',xTicks); set(hAx(2),'XTick',[]);
    set(hAx(1),'XTickMode','Manual'), set(hAx(2),'XTickMode','Manual')
    set(hAx(1),'YTickMode','Manual'), set(hAx(2),'YTickMode','Manual')
  end
  xlabel(unit,'FontSize',fs)
%
% Reset the default text interpreter.
%
  set(0,'DefaultTextInterpreter',original_interpreter);
%
end

% Copyright 2016 Kevin J. Kircher. See bldg/license.txt for full license
% information.