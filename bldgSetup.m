function bldgSetup
%BLDGSETUP sets up and tests the BLDG distribution.

%
% Get the directory delimiter.
%
  wd = pwd;
  delimiter = wd(end-4);
  last_word = wd(end-3:end);
  if ~strcmp(last_word,'bldg')
    error('Please run bldgSetup from the bldg directory.')
  end
%
% Try to add the required folders to the MATLAB path.
%
  fprintf('\nAdding required folder to your path... ')
  addFailed = 0;
  try
    addpath(strcat(pwd,delimiter,'functions'))
  catch
    addFailed = 1;
  end
%
% If that failed, print alternate instructions and return.
%
  if addFailed
    fprintf('\nError: bldgSetup failed, likely due to permissions.\n')
    fprintf('Please manually add the subfolder `functions` to your path.\n')
    fprintf('Once this is complete, save your path and BLDG should be ready to go.\n')
    return
  else
    fprintf('success.\n')
  end
%
% Try to save the MATLAB path.
%
  fprintf('Saving your path... ')
  saveFailed = savepath;
  if saveFailed
    fprintf('\nWarning: bldgSetup failed to save your path, likely due to permissions.\n')
    fprintf('Please save your path manually.\n')
    fprintf('Otherwise, bldgSetup will need to be run every time MATLAB restarts.\n')
  else
    fprintf('success.\n')
  end
%
% Define the timing and input signals for a step change in the outdoor air
% temperature, with no other forcing.
%
  nt = 1e3; t = linspace(0,24*3600,nt)';
  W = zeros(nt,6);
  W(:,1) = 273*ones(nt,1);
  U = zeros(nt,1);
  N = 51;
  x0 = 293*ones(N+2,1);
%
% Try to define a default bldg object and simulate it.
%
  fprintf('Trying a simple example... ')
  try
    b = bldg;
    X = bsim(b,t,W,x0,U);
  catch
    fprintf('Error: the installation failed for unknown reasons.\n')  
    fprintf('Please see the Support section of the BLDG Users` Guide to get help.')  
  end
  fprintf('success.\nBLDG is ready to go!\n\n')
%
end

