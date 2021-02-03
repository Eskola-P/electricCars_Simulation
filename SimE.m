
% ------------------------------------- %
% Generate Simulations for              % 
% - Global Irradiation on the ground    %
% - Global Irradiation on tilted ground %
% - Total Irradiation on tilted ground  % 
%                                       %
% Version 2.0                           %
%                                       %
% For scientific use only               %
% Copyright: Stephan Schlüter           %
% Last update: 31.01.2021               %
% ------------------------------------- %


function E_tot = simE(simsCloud,simsTemp,simsTimeStamp,data_geo)

%year = 2020
%simsCloud = cloudSim
%simsTemp = tempSims
%simsTimeStamp = TimeStamp 

  
% Simulate E^glo,tilt 30% southwards

% parameters for estimating E^glo,tilt depending on the hour of day

% hour intercept	c_t	(90 - theta)	E^Extra	T^air

pars_tilt =  [0.3930 	 0.0063 	-0.2051 	 0.0089 	-0.0011 	-0.1730;
              0.8892	-0.0336 	-0.0559 	 0.0017 	 0.0007  	 0.0954;
              0.9814 	-0.0427 	-0.0437 	 0.0013 	 0.0044 	 0.1819;
              1.0066 	-0.0452 	-0.0423 	 0.0014 	 0.0060 	 0.2181;
              0.7323 	-0.0451 	-0.0583 	 0.0027 	 0.0060 	 0.1548; 
              0.7920 	-0.0403 	-0.0414 	 0.0018 	 0.0117 	 0.1808; 
              0.9114 	-0.0354 	-0.0332 	 0.0012 	 0.0130 	 0.2369; 
              0.9308 	-0.0398 	-0.0267 	 0.0009 	 0.0119 	 0.1983;
              0.8582 	-0.0449 	-0.0334 	 0.0013 	 0.0106 	 0.1491; 
              0.6622 	-0.0401 	-0.0485 	 0.0021 	 0.0078 	 0.1243; 
              0.2767 	-0.0271 	-0.1058 	 0.0050 	 0.0045 	 0.0870; 
              0.1091 	-0.0173 	-0.1230 	 0.0058 	 0.0018 	 0.0498; 
              0.0160 	-0.0060 	-0.0104 	 0.0009 	 0.0005 	 0.0273;
             -0.0280  0.0007 	-0.0097 	 0.0009 	-0.0007 	 0.0174];


E_glo_tilt = zeros(size(simsCloud));
E_glo      = zeros(size(simsCloud));
E_refl     = zeros(size(simsCloud));

for k = 1:size(simsCloud,1)
    
    % Epsilon ############
    if (simsTimeStamp(k,4) <6 | simsTimeStamp(k,4) > 19 )
        pars_k = [0 0 0 0 0 0];
    else
        % we just start at 6 a.m., that's why -5
        pars_k = pars_tilt(simsTimeStamp(k,4)-5,:);
    end
    
    factor1 = pars_k(1);
    factor2 = pars_k(2)*simsCloud(k,:);
    factor3 = pars_k(3)*(90 - (acos(data_geo(k,2))*180/pi));
    factor4 = pars_k(4)*data_geo(k,6);
    factor5 = pars_k(5)*simsTemp(k,2:end);
    factor6 = pars_k(6)*(-cos(simsTimeStamp(k,2)*2*pi/12));
    epsilon = factor1 + factor2 + factor3 + factor4 + factor5 + factor6;
    
    E_glo_tilt(k,:) = data_geo(k,6).*epsilon;
    
end


% ################################# #
% E^tot,tilt = E^glo,tilt + E^refl
% ################################# #


% hour intercept	c_t	(90 - theta)	E^Extra	T^air
pars_glo =  [0.8664 	-0.0447 	-0.1085 	 0.0041 	 0.0052 	 0.0565 ;
             0.9336     -0.0466 	-0.1556 	 0.0061	     0.0038	     0.1519 ;
             1.1537 	-0.0460 	-0.0251 	 0.0002	     0.0044 	 0.2764 ;
             1.1165 	-0.0449 	-0.0077 	-0.0003	     0.0061      0.2492 ;
             0.8488 	-0.0423 	-0.0356 	 0.0014	     0.0060      0.1797 ;
             0.7877 	-0.0377 	-0.0230 	 0.0009	     0.0083      0.1552 ;
             0.7460 	-0.0325 	-0.0197 	 0.0008	     0.0082      0.1694 ;
             0.7132 	-0.0306 	-0.0154 	 0.0006	     0.0078      0.1571 ;
             0.6043     -0.0285 	-0.0189 	 0.0008      0.0071      0.1300 ;
             0.5085 	-0.0232 	-0.0219 	 0.0009      0.0062      0.1297 ;
             0.2547 	-0.0184 	-0.0654 	 0.0031      0.0042 	 0.1147 ;
             0.0964 	-0.0134 	-0.0141 	 0.0010      0.0020      0.0379 ;
             0.0256 	-0.0084 	-0.0180 	 0.0013 	 0.0008 	 0.0592 ;
            -0.0276 	-0.0044 	-0.0201 	 0.0015 	-0.0005 	 0.0320 ];





for k = 1:size(simsCloud,1)

  % Epsilon #############
  if (simsTimeStamp(k,4) <6 | simsTimeStamp(k,4) > 19)   
    % before 6 and after 7 p.m. we have no irradiation
    pars_k = [0 0 0 0 0 0];
  else  
    % we just start at 6 a.m., that's why -5
    pars_k = pars_glo(simsTimeStamp(k,4)-5,:);
  end
 
  E_extra = data_geo(k,6).*ones(1,size(simsCloud,2));
  
  factor1 = pars_k(1);
  factor2 = pars_k(2)*simsCloud(k,:);
  factor3 = pars_k(3)*(90 - (acos(data_geo(k,2))*180/pi));
  factor4 = pars_k(4)*E_extra;
  factor5 = pars_k(5)*simsTemp(k,2:end);
  factor6 = pars_k(6)*(-cos(simsTimeStamp(k,2)*2*pi/12));
  epsilon = factor1 + factor2 + factor3 + factor4 + factor5 + factor6;
  
  E_glo(k,:) = data_geo(k,6).*epsilon;
  end


% compute the Albedo for the reflexive irradiation
% angle = theta_z
angle = (90 - (acos(data_geo(:,2)).*180/pi));
albedo = (angle > 30).*20 + (angle <= 30).*(100 - (8/3).*angle);
coeff = (albedo./200).*(1 - cos(30*pi/180));

for k = 1:size(E_glo,1)
  E_refl(k,:) = E_glo(k,:)*coeff(k);
end


E_tot = E_glo_tilt + E_refl;

% correct: if theta_z < 0 - everything is zero 
% alternatively, if cos(theta_z)/cos(theta) < 0 --> sun does not shine onto the panel
for k = 1:size(E_glo,1)
  if (angle(k) <= 0 | cos(data_geo(k,5))/data_geo(k,2) < 0) 
    E_tot(k,:) = zeros(1,size(E_tot,2));
  end
end

