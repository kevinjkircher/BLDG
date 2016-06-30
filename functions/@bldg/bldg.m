classdef bldg
%BLDG is a class that represents the simplest building by its geometry and 
%material properties. SI units are used throughout.
  
  properties
  %%% Geometry %%%     
%     l = 0.6;  							          % wall thickness
    l = 0.1;  							          % wall thickness
    lg = 0.0032;						          % glass thickness
%     la = 1;                           % room length
    la = 7;                           % room length
%     A                                 % wall/glass surface area
    A = 21;                           % wall/glass surface area
    gam = pi;                         % wall azimuth angle in radians
    phi = 0.7118;						          % latitude in radians
    
  %%% Materials %%%
%     a = 5e-8;                         % wall thermal diffusivity
    a = 7.5e-8;                       % wall thermal diffusivity
    k = 0.4;                          % wall thermal conductivity
    eps = 0.9;                        % wall longwave emissivity
%     as = 0.4;                         % wall shortwave absorptivity
    as = 0.2;                         % wall shortwave absorptivity
    Cg                                % glass thermal capacitance
    epsg = 0.2;                       % glass longwave emissivity
    ng = 1.53;                        % glass index of refraction
    mug = 19.6;                       % glass attenuation coefficient
    Ca                                % air thermal capacitance
        
  %%% Internal heat sources %%%
    mdot                              % infiltration mass flow rate
    zc = 1;                           % control systems convective fraction
    zp = 0.4;                         % people convective fraction
    zl = 0.75;                        % lights convective fraction
    ze = 0.2;                         % equipment convective fraction
    eta = 0.2;							          % lighting efficiency
  end
  
  properties (Hidden, Constant)
    sig = 5.67e-8;                    % Stefan-Boltzmann constant
    rhoa = 1.225;                     % air density
    ca = 1005;                        % air specific heat
  end
  
  properties (Hidden)
%     areaRatio = 10;                   % sqrt(A)/la
    agBar = 0.0682;                   % glass diffuse absorptivity
    taugBar = 0.7786;                 % glass diffuse transmissivity
  end
  
  methods
  %%% Constructor %%%
    function b = bldg()
    %This function is the constructor for the bldg class.
          
    %%% Geometry %%%
%       b.A = (b.areaRatio*b.la)^2;         % wall/glass area
      
    %%% Materials %%%
      rhog = 2530;                        % glass density
      cg = 880;                           % glass specific heat
      b.Cg = rhog*cg*b.A*b.lg;            % glass thermal capacitance
%       capacitanceRatio = 20;              % Ca/(rhoa*ca*A*la)
      capacitanceRatio = 10;              % Ca/(rhoa*ca*A*la)
      b.Ca = capacitanceRatio*...
        b.rhoa*b.ca*b.A*b.la;             % air thermal capacitance
    
    %%% Internal heat sources %%%
      ACH = 2;
      b.mdot = ACH*b.rhoa*b.la*b.A/3600;  % infiltration mass flow rate
    end
  
  %%% Geometry setters %%%
    function b = set.l(b,l)
    %This function sets the wall thickness.
      if l <= 0
        error('Wall thickness must be positive.')
      else
        b.l = l;
      end
    end
    
    function b = set.lg(b,lg)
    %This function sets the glass thickness.
      if lg <= 0
        error('Glass thickness must be positive.')
      else
        b.lg = lg;
        b = updateAbsTrans(b);
      end
    end
    
    function b = set.la(b,la)
    %This function sets the room length.
      if la <= 0
        error('Room length must be positive.')
      else
        b.la = la;
%         b = updateAreaRatio(b);
      end
    end
    
    function b = set.A(b,A)
    %This function sets the wall/glass area.
      if A <= 0
        error('Wall/glass area must be positive.')
      else
        b.A = A;
%         b = updateAreaRatio(b);
      end
    end
    
    function b = set.gam(b,gam)
    %This function sets the wall azimuth angle.
      if gam < 0 || gam >= 2*pi
        error('Wall azimuth angle must be in [0, 2*pi) radians.')
      else
        b.gam = gam;
      end
    end
    
    function b = set.phi(b,phi)
    %This function sets the building latitude.
      if phi < -pi/2 || phi > pi/2
        error('Building latitude must be in [-pi/2, pi/2] radians.')
      else
        b.phi = phi;
      end
    end
    
  %%% Material setters %%%
    function b = set.a(b,a)
    %This function sets the wall thermal diffusivity.
      if a <= 0
        error('Wall thermal diffusivity must be positive.')
      else
        b.a = a;
      end
    end
    
    function b = set.k(b,k)
    %This function sets the wall thermal conductivity.
      if k <= 0
        error('Wall thermal conductivity must be positive.')
      else
        b.k = k;
      end
    end
    
    function b = set.eps(b,eps)
    %This function sets the wall longwave emissivity.
      if eps <= 0 || eps > 1
        error('Wall longwave emissivity must be in (0,1].')
      else
        b.eps = eps;
      end
    end
    
    function b = set.as(b,as)
    %This function sets the wall shortwave absorptivity.
      if as <= 0 || as > 1
        error('Wall shortwave absorptivity must be in (0,1].')
      else
        b.as = as;
      end
    end
    
    function b = set.Cg(b,Cg)
    %This function sets the glass thermal capacitance.
      if Cg <= 0
        error('Glass thermal capacitance must be positive.')
      else
        b.Cg = Cg;
      end
    end
    
    function b = set.epsg(b,epsg)
    %This function sets the glass longwave emissivity.
      if epsg <= 0 || epsg > 1
        error('Glass longwave emissivity must be in (0,1].')
      else
        b.epsg = epsg;
      end
    end
    
    function b = set.ng(b,ng)
    %This function sets the glass index of refraction.
      if ng <= 0
        error('Glass index of refraction must be positive.')
      else
        b.ng = ng;
        b = updateAbsTrans(b);
      end
    end
    
    function b = set.mug(b,mug)
    %This function sets the glass attenuation coefficient.
      if mug <= 0
        error('Glass attenuation coefficient must be positive.')
      else
        b.mug = mug;
        b = updateAbsTrans(b);
      end
    end
    
    function b = set.Ca(b,Ca)
    %This function sets the room air thermal capacitance.
      if Ca <= 0
        error('Air thermal capacitance must be positive.')
      else
        b.Ca = Ca;
      end
    end
    
  %%% Internal heat source setters %%%
    function b = set.mdot(b,mdot)
    %This function sets the outdoor air infiltration rate.
      if mdot < 0
        error('Outdoor air infiltration rate must be nonnegative.')
      else
        b.mdot = mdot;
      end
    end
    
    function b = set.zc(b,zc)
    %This function sets the convective fraction of heat from control 
    %systems.
      if zc <= 0 || zc > 1
        error('Convective fraction of heat from control systems must be in (0,1].')
      else
        b.zc = zc;
      end
    end
    
    function b = set.zp(b,zp)
    %This function sets the convective fraction of occupant body heat.
      if zp < 0 || zp > 1
        error('Convective fraction of occupant body heat must be in [0,1].')
      else
        b.zp = zp;
      end
    end
    
    function b = set.zl(b,zl)
    %This function sets the convective fraction of heat from lights.
      if zl < 0 || zl > 1
        error('Convective fraction of heat from lights must be in [0,1].')
      else
        b.zl = zl;
      end
    end
    
    function b = set.ze(b,ze)
    %This function sets the convective fraction of heat from lights.
      if ze < 0 || ze > 1
        error('Convective fraction of heat from miscellaneous equipment must be in [0,1].')
      else
        b.ze = ze;
      end
    end
    
    function b = set.eta(b,eta)
    %This function sets the convective fraction of heat from lights.
      if eta <= 0 || eta >= 1
        error('Lighting efficiency must be in (0,1).')
      else
        b.eta = eta;
      end
    end
    
  %%% Hidden property updaters %%%
%     function b = updateAreaRatio(b)
%     %This function updates the ratio of the wall/glass height (or width) to
%     %the room length.
%       b.areaRatio = sqrt(b.A)/b.la;
%       if b.areaRatio < 8
%         warning_state = warning;          % get original warning state
%         warning off backtrace
%         warning('Model assumptions are inaccurate when sqrt(A)/la is much less than 10.')
%         warning(warning_state)            % return to original warning state
%       end
%     end
    
    function b = updateAbsTrans(b)
    %This function updates the glass diffuse absorptivity and
    %transmissivity.
      nthg = 1e3; thgs = linspace(0,pi/2,nthg);
      [ags,taugs] = getAbsTrans(thgs,b.lg,b.ng,b.mug);
      b.agBar = 2*trapz(thgs, ags.*sin(thgs).*cos(thgs));
      b.taugBar = 2*trapz(thgs, taugs.*sin(thgs).*cos(thgs));
    end
    
  %%% Display %%%
    function disp(b)
    %This function displays the contents of a bldg object.
    %%% Geometry %%%
      fprintf('Object of class BLDG.\n')
      fprintf('\nGeometric parameters:\n')
      fprintf('  Wall thickness: l = %.3g m.\n',b.l)
      fprintf('  Glass thickness: lg = %.3g m.\n',b.lg)
      fprintf('  Room length: la = %.3g m.\n',b.la)
      fprintf('  Wall/glass surface area: A = %.3g m^2.\n',b.A)
      fprintf('  Wall azimuth (0 = north): gam = %.3g rad.\n',b.gam)
      fprintf('  Building latitude (pi/2 = north pole): phi = %.3g rad.\n',b.phi)
      
    %%% Materials %%%
      fprintf('\nMaterial properties:\n')
      fprintf('  Wall diffusivity: a = %.3g m^2/s.\n',b.a)
      fprintf('  Wall conductivity: k = %.3g W/(m*K).\n',b.k)
      fprintf('  Wall longwave emissivity: eps = %.3g.\n',b.eps)
      fprintf('  Wall shortwave absorptivity: as = %.3g.\n',b.as)
      fprintf('  Glass capacitance: Cg = %.3g J/K.\n',b.Cg)
      fprintf('  Glass longwave emissivity: epsg = %.3g.\n',b.epsg)
      fprintf('  Glass index of refraction: ng = %.3g.\n',b.ng)
      fprintf('  Glass attenuation coefficient: mug = %.3g m^(-1).\n',b.mug)
      fprintf('  Room capacitance: Ca = %.3g J/K.\n',b.Ca)
      
    %%% Internal heat sources %%%
      fprintf('\nInternal heat source properties:\n')
      fprintf('  Outdoor air infiltration mass flow rate: mdot = %.3g kg/s.\n',b.mdot)
      fprintf('  Convective fraction of heat from control systems: zc = %.3g.\n',b.zc)
      fprintf('  Convective fraction of occupant body heat: zp = %.3g.\n',b.zp)
      fprintf('  Convective fraction of heat from lighting: zl = %.3g.\n',b.zl)
      fprintf('  Convective fraction of heat from miscellaneous equipment: ze = %.3g.\n',b.ze)
      fprintf('  Lighting efficiency: eta = %.3g.\n',b.eta)
    end
    
  %%% Methods in separate files %%%%
    [b,weather] = importWeather(b,fileName,dt)
    gains = generateGains(b,t,varargin)
    [X,Q] = bsim(b,t,W,x0,varargin)
    [qgm,qgp,qm,qp,qam,qap,Qa] = getFluxes(b,x,Tinf,Id,Ibp,Ibg,thg,Qp,Ql,Qe,Qc)
    Qc = getLoads(b,X,W,t)
    [Tinf,Id,Ibp,Ibg,thg] = getWeather(b,W,t,tNum)
    x0 = precondition(b,t,weather,gains,N,Ts)
  end
end

% Copyright 2016 Kevin J. Kircher. See bldg/license.txt for full license
% information.