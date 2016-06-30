%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Examples from the 'Quick start' section of the BLDG Users' Guide. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.1 Defining a building %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define a default bldg object and display its properties.
b = bldg

% Redefine the wall shortwave absorptivity.
b.as = 0.1;

% Try to set an unphysical wall length. Display the error.
try
  b.l = -1;
catch err
  fprintf('When asked to make the wall thickness negative, BLDG threw the following error: \n\t')
  disp(err.message)
end

try input('Press Enter/Return for the next example...'), clc, catch, end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.2 Step response to outdoor temperature %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define a 5-day simulation time span with 15-minute time steps.
dt = 15*60;                                       % time step in seconds
t = 0:dt:5*24*3600;

% Define the disturbance and control signals.
nt = length(t); 
nw = 6;                                           % disturbance dimension
nu = 1;                                           % control dimension
W = zeros(nt,nw);
W(:,1) = 273;                                     % outdoor air temperature
W(1:round(nt/10),1) = 293;
U = zeros(nt,nu);

% Define the number of wall nodes and the initial state.
N = 50;
x0 = 293*ones(N+2,1);

% Simulate the step response and plot the results.
X = bsim(b,t,W,x0,U);
figure(1), plotResults(t,X,W,U)
%dir = '/Users/Kevin/Documents/Research/Projects/Simplest building/User guide/Graphics/'
%print('-depsc2', strcat(dir, 'oat_step.eps'))

try input('Press Enter/Return for the next example...'), clc, catch, end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.2 (continued) Step response to control heat flow %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Update the disturbance and control signals and the initial state.
W(:,1) = 273*ones(nt,1);
Qc0 = 5e3;
U(round(nt/10):end) = Qc0;
x0 = 273*ones(N+2,1);

% Simulate the step response and plot the results.
X = bsim(b,t,W,x0,U);
figure(1), plotResults(t,X,W,U)
%print('-depsc2', strcat(dir, 'control_step.eps'))

try input('Press Enter/Return for the next example...'), clc, catch, end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.3 Simulating a forced-air heating system %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the supply air temperature.
Tc = 325; 

% Define the parameters of the sinusoidal supply air mass flow rate.
mdotc_mean = 0.1;
mdotc_amp = mdotc_mean; 
mdotc_period = 4*3600;

% Update the control signal to include supply air mass flow rates in column
% 1 and supply air temperatures in column 2.
nu = 2;
U = zeros(nt,nu);
U(:,2) = Tc;
U(:,1) = mdotc_mean + mdotc_amp*sin((2*pi/mdotc_period)*t)';

% Simulate the transient response and plot the results.
X = bsim(b,t,W,x0,U);
figure(1), plotResults(t,X,W,U)
%print('-depsc2', strcat(dir, 'forced_air.eps'))

try input('Press Enter/Return for the next example...'), clc, catch, end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.4 Calculating peak heating and cooling loads %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Redefine the time step to one hour.
dt = 3600;

% Import a year of weather from a TMY3 file for New York City.
[b,weather] = importWeather(b,'NYC_TMY3.csv',dt);

% Generate a year of internal heat flows.
gains = generateGains(b,weather.tw,'residential');

% Set the simulation time span to the full year for which the weather data
% is available.
t = weather.tw;

% Find an initial state consistent with the weather and internal gains.
Ts = 293;                                         % temperature setpoint
x0 = precondition(b,t,weather,gains,N,Ts);

% Pack the disturbances into the matrix form accepted by bsim.
W = getDisturbances(weather,gains,t);

% Compute the heating and cooling loads by simulating the building under
% perfect control.
[X,U] = bsim(b,t,W,x0);
U_max = max(U);
U_min = min(U);

% Display the maximum heating and cooling loads.
fprintf('Peak heating load: %.3g kW.\n',U_max/1e3)
fprintf('Peak cooling load: %.3g kW.\n',max(0,-min(U)/1e3))

try input('Press Enter/Return for the next example...'), clc, catch, end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.5 Simulating closed-loop operation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the heater capacity.
Qc0 = 2*U_max;

% Shorten the time step to 10 minutes and update the weather and internal
% gains.
dt = 5*60;
[b,weather] = importWeather(b,'NYC_TMY3.csv',dt);
gains = generateGains(b,weather.tw,'residential');

% Define a 1-day simulation time span starting at midnight on January 18.
nd0 = 18;
t0 = 0; tf = 24*3600;
t = getTiming(weather.tw,nd0,t0,tf);

% Get the disturbances and initial state.
W = getDisturbances(weather,gains,t);
x0 = precondition(b,t,weather,gains,N,Ts);

% Initialize the state and control matrices.
nt = length(t);
X = zeros(nt,N+2);
X(1,:) = x0';
nu = 1;
U = zeros(nt,nu);

% Define the thermostat deadband temperatures.
T_high = Ts + 2;
T_low = Ts - 2;

% Loop through the simulation time span.
tic
for k = 1 : nt-1

  % Let the system evolve, with bsim acting as a nonlinear, time-varying, 
  % discrete-time dynamics function.
  x_new = bsim(b,t(k:k+1),W(k,:),X(k,:),U(k));

  % Decide the control.
  Ta = x_new(end);
  if Ta < T_low
    U(k+1) = Qc0;
  elseif Ta > T_high
    U(k+1) = 0;
  else
    U(k+1) = U(k);
  end

  % Store the data from this iteration.
  X(k+1,:) = x_new';

end
closedLoopTime = toc;

% Plot the simulation inputs and outputs.
figure(1), plotInputs(t,W)
%print('-depsc2', strcat(dir, 'closed_loop_inputs.eps'))
figure(2), plotResults(t,X,W,U,T_low,T_high)
%print('-depsc2', strcat(dir, 'closed_loop_results.eps'))

% Compare the open-loop and closed-loop simulation times.
tic
X = bsim(b,t,W,x0,Qc0/2*ones(nt,1));
openLoopTime = toc;
fprintf('Closed-loop simulation was %.3g times slower than open-loop.\n',...
  closedLoopTime/openLoopTime)

