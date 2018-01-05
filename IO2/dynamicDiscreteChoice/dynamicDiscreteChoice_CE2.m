rng(1)

nPeriods = 100
nFirms = 1000

tolFixedPoint = 1e-10
			
nSuppX = 5;
supportX = (1:nSuppX)'
capPi = 1./(1+abs(ones(nSuppX,1)*(1:nSuppX)-(1:nSuppX)'*ones(1,nSuppX)));
capPi = capPi./(sum(capPi')'*ones(1,nSuppX))
beta = [-0.1*nSuppX;0.2]
delta = [0;1]
rho = 0.95	

[u0,u1] = flowpayoffs(supportX,beta,delta); 
[capU0,capU1] = fixedPoint(u0,u1,capPi,rho,tolFixedPoint,@bellman,[],[]);
deltaU = capU1-capU0;

[choices,iX] = simulateData(deltaU,capPi,nPeriods,nFirms);

% tolFixedPoint = 1e-10                  
% tolFixedPoint = 1e-1                                  % w: change inner loop tol here !
tolFixedPoint = 10                                      % w: when the inner loop tol is large the estimated results change

objectiveFunction = @(parameters)negLogLik(choices,iX,supportX,capPi,parameters(1:2),[delta(1);parameters(3)],...
                                           rho,@flowpayoffs,@bellman,@fixedPoint_CE2,tolFixedPoint)

startvalues = [-1;-0.1;0.5];

lowerBounds = -Inf*ones(size(startvalues));
lowerBounds(3) = 0;

% w: change outer loop tol here
% OptimizerOptions = optimset('Display','iter','Algorithm','interior-point','AlwaysHonorConstraints','bounds',...
%                             'GradObj','on','TolFun',1E-6,'TolX',1E-10,'DerivativeCheck','off','TypicalX',[beta;delta(2)]);
OptimizerOptions = optimset('Display','iter','Algorithm','interior-point','AlwaysHonorConstraints','bounds',...
                            'GradObj','on','TolFun',1E-1,'TolX',1E-1,'DerivativeCheck','off','TypicalX',[beta;delta(2)]);
% OptimizerOptions = optimset('Display','iter','Algorithm','interior-point','AlwaysHonorConstraints','bounds',...
%                             'GradObj','on','TolFun',1E-14,'TolX',1E-14,'DerivativeCheck','off','TypicalX',[beta;delta(2)]);
[maxLikEstimates,~,exitflag] = fmincon(objectiveFunction,startvalues,[],[],[],[],lowerBounds,[],[],OptimizerOptions)

[~,~,informationMatrix] = objectiveFunction(maxLikEstimates);
standardErrors = diag(sqrt(inv(informationMatrix)));

disp('Summary of Results');
disp('--------------------------------------------');
disp('      true     start     estim      ste.');
disp([[beta;delta(2)] startvalues maxLikEstimates standardErrors]);

piHat = estimatePi(iX,nSuppX)
