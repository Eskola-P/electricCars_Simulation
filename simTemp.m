% ------------------------------------- %
% Generate Temperature Simulations      % 
% Version 2.0                           %
%                                       %
% For scientific use only               %
% Copyright: Stephan Schlüter           %
% Last update: 31.01.2021               %
% ------------------------------------- %

function erg = simTemp(X,beginDate,endDate,nSims)

% X: [year month day hour temperature]
% we expect hourly data

% -------------------------------------------------------------------------
% Model Calibration
a = X(:,end);

stabw = std(a);
a = a./std(a);

idx = (1:size(a,1))';

%Create function
fun = @(x,xdata)(x(1) + x(2).*sin(xdata.*(2*pi)./(365.25*24) - x(3).*(2*pi)./(365.25*24) ));
%Create an initial guess.
x0 = [0,0.8,0];
%Solve the bounded fitting problem.
%x = lsqcurvefit(fun,x0,xdata,X(:,4),lb,ub)
pars = lsqcurvefit(fun,x0,idx,a);

trend = (pars(1) + pars(2).*sin(idx.*(2*pi)./(365.25*24) - pars(3).*(2*pi)./(365.25*24) )).*stabw;


vec1 = trend(24:(size(X,1)-1)) - X(24:(size(X,1)-1),end);
vec2 = trend(23:(size(X,1)-2)) - X(23:(size(X,1)-2),end);
vec3 = trend(1:(size(X,1)-24)) - X(1:(size(X,1)-24),end);

Z = [ vec1 vec2 vec3];
y = X(25:(size(X,1)),end) - X(24:(size(X,1)-1),end);

pars_LSQ = inv(Z'*Z)*Z'*y;

residuals = y - pars_LSQ(1).*vec1 - pars_LSQ(2).*vec2 - pars_LSQ(3).*vec3;

stabwTemp = std(residuals);


% -------------------------------------------------------------------------
% Model Simulation

% Trend ausrollen

% t1 = datetime(year(x(1,1)),month(x(1,1)),day(x(1,1)),hour(x(1,1)),0,0);
t1 = datetime(X(1,1),X(1,2),X(1,3),X(1,4),0,0);
t2 = datetime(year(beginDate),month(beginDate),day(beginDate),hour(beginDate),0,0);
t3 = datetime(year(endDate),month(endDate),day(endDate),hour(endDate),0,0);

indexVec = (datenum(t1:hours(1):t3))';
indexVec(:,2) = 1:length(indexVec);
trend_all = (pars(1) + pars(2).*sin(indexVec(:,2).*(2*pi)./(365.25*24) - pars(3).*(2*pi)./(365.25*24) )).*stabw;

% Simulation ab letzten Zeitpunkt von X, muss dann entsprechend den
% Simulationszeitraum final ausschneiden.

% Achtung, weil ich 24h-lag drin habe, brauche ich noch die 24 letzten
% werte
sims = zeros(length(trend_all),nSims+1);
sims(:,1) = indexVec(:,1);
sims(1:size(X,1),2:end) = repmat(X(:,end),1,nSims);

for k = (size(X,1)+1):size(sims,1)
     zufall = normrnd(0,stabwTemp,1,nSims);
     sims(k,2:end) = sims(k-1,2:end) + pars_LSQ(1).*(trend_all(k-1) - sims(k-1,2:end)) ...
         + pars_LSQ(2).*(trend_all(k-2) - sims(k-2,2:end)) ...
         + pars_LSQ(3).*(trend_all(k-24) - sims(k-24,2:end)) ...
         + zufall;
end


% Muss nur noch den entsprechenden Zeitraum ausschneiden.

anfang = find(sims(:,1) == beginDate);
ende   = find(sims(:,1) == endDate);

erg = sims(anfang:ende,:);
