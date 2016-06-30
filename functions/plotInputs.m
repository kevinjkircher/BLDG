function plotInputs(t,W)
%PLOTINPUTS plots the exogenous inputs to the simplest building.

% Inputs:
%	t, an array of nt simulation times in solar seconds.
%	W, a nt x 5 matrix containing the following disturbance histories:
%		W(:,1) = Tinf, the outdoor air temperature.
%		W(:,2) = Ih, the total horizontal irradiance.
%		W(:,3) = Ibn, the beam normal irradiance.
%		W(:,4) = Qp, theinternal heat gained from people.
%		W(:,5) = Ql, the internal heat gained from lights.
%		W(:,6) = Qe, the internal heat gained from other equipment.

%
% Get the timing.
%
  nd0 = floor(t(1)/(24*3600)) + 1;					% initial day number
  s = t - (nd0-1)*24*3600;							% seconds since midnight
  h = s/3600;										% hours since midnight
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
% Get the inputs.
%
  Tinf = W(:,1); Ih = W(:,2); Ibn = W(:,3);
  [Qp,Ql,Qe] = getGains(W,t,t);
%
% Get the default text interpreter. Temporarily change it to LateX.
%
  original_interpreter = get(0,'DefaultTextInterpreter');
  set(0,'DefaultTextInterpreter','LaTeX');
%
% Plot the weather in the first subplot.
%
  clf, subplot(2,1,1)
  [hAx,h1,h2] = plotyy(x,Tinf - 273, x,[Ih,Ibn]);
  set(h1,'linewidth',2); set(h2,'linewidth',2);
  set(h2(1),'color','k'); set(h2(2),'color','r')
  leg = legend('$T_\infty$','$I_h$','$I_\perp^b$',...
    'Location','NorthOutside','Orientation','Horizontal');
  set(leg,'color','w','FontSize',fs+2,'Interpreter','LaTeX')
  ylabel(hAx(1),'Temperature (C)','FontSize',fs)
  ylabel(hAx(2),'Irradiance (W/m$^2$)','FontSize',fs)
  xlim(hAx(1),xLimits), xlim(hAx(2),xLimits)
  yLimits2 = get(hAx(2),'YLim'); yLimits2(2) = yLimits2(1) + 1.2*(yLimits2(2) - yLimits2(1));
  ylim(hAx(2),yLimits2)
  set(hAx(1),'XTick',xTicks); set(hAx(2),'XTick',xTicks);
  set(hAx(1),'XTickMode','Manual','YTickMode','Manual')
  set(hAx(2),'XTickMode','Manual','YTickMode','Manual')
%
% Plot the internal gains in the second subplot.
%
  subplot(2,1,2), hold on
  plot(x,Qp,'k',x,Ql,'m', x,Qe,'g','linewidth',2)
  yLimits = get(gca,'YLim'); yLimits(2) = yLimits(1) + 1.1*(yLimits(2) - yLimits(1));
  ylim(yLimits), box on
  leg = legend('$Q_p$','$Q_l$','$Q_e$',...
    'Location','NorthOutside','Orientation','Horizontal');
  set(leg,'FontSize',fs+2,'Interpreter','LaTeX')
  xlabel(unit,'FontSize',fs), xlim(xLimits)
  ylabel('Internal gains (W)','FontSize',fs)
  set(gca,'XTick',xTicks)
  set(gca,'XTickMode','Manual'), set(gca,'YTickMode','Manual')
%
% Reset the default text interpreter.
%
  set(0,'DefaultTextInterpreter',original_interpreter);
%
end

% Copyright 2016 Kevin J. Kircher. See bldg/license.txt for full license
% information.