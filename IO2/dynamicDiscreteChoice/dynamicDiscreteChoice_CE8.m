clc
clear

%% NFXP
rng(1)                         % w: set seed for Monte Carlo experiment
%{
\section{The Script that Puts It All Together\label{script}}

	The script in |dynamicDiscreteChoice.m| simulates data and computed maximum likelihood estimates using the nested fixed point (NFXP) method of \cite{ecta87:rust} and \cite{nh94:rust}. It takes ${\cal
	X},\delta_0,\rho$ to be known, either takes $\Pi$ to be known or estimates it in a first stage, and focuses on maximum partial likelihood estimation of the remaining parameters; $\beta_0,\beta_1,\delta_1$; from conditional choice probabilities. 
	
\subsection{Simulating Data}
	
	First, we set the number of time periods (|nPeriods|) and firms (|nFirms|) that we would like to have in our sample.
%}
nPeriods = 100
nFirms = 1000
%{
	We also set the tolerance |tolFixedPoint| on the fixed point $U$ of $\Psi$ that we will use to determine the simulation's entry and exit rules. This same tolerance will also be used when solving the model in the inner loop of the NFXP procedure.
%}
tolFixedPoint = 1e-10
%{
Next, we specify the values of the model's parameters used in the simulation: 
	\begin{dictionary}
	\item{|nSuppX|} the scalar number $K$ of elements of ${\cal X}$;
	\item{|supportX|} the $K\times 1$ vector ${\cal X}$ with the support points of $X_t$;	
	\item{|capPi|} the $K\times K$ Markov transition matrix $\Pi$ for $\{X_t\}$, with typical element $\Pi_{ij}=\Pr(X_{t+1}=x^j|X_t=x^i)$;
	\item{|beta|} the $2\times 1$ vector $\beta$ with the parameters of the flow profit of active firms;
	\item{|delta|} the $2\times 1$ vector of exit and entry costs $\delta$; and
	\item{|rho|} the scalar discount factor $\rho$.
	\end{dictionary}
	%}													
nSuppX = 5;
supportX = (1:nSuppX)'
capPi = 1./(1+abs(ones(nSuppX,1)*(1:nSuppX)-(1:nSuppX)'*ones(1,nSuppX)));
capPi = capPi./(sum(capPi')'*ones(1,nSuppX))
beta = [-0.1*nSuppX;0.2]
delta = [0;1]
rho = 0.95	
%{
For these parameter values, we compute the flow payoffs $u_0$ (|u0|) and $u_1$ (|u1|), the choice-specific expected discounted values $U_0$ (|capU0|) and $U_1$ (|capU1|), and their contrast $\Delta U$ (|deltaU|).
%}
[u0,u1] = flowpayoffs(supportX,beta,delta); 
[capU0,capU1] = fixedPoint(u0,u1,capPi,rho,tolFixedPoint,@bellman,[],[]);
deltaU = capU1-capU0;
%{
	With $\Delta U$ computed, and $\Pi$ specified, we proceed to simulate a $T\times N$ matrix of choices |choices| and a $T\times N$ matrix of states |iX| (recall from Section \ref{simulate} that |iX| contains indices that point to elements of ${\cal X}$ rather than those values themselves).
%}

nSimu = 10;                                % w: set # of simulated data sets for running the Monte Carlo experiment

maxLikEstimates_mat = zeros(3,nSimu);
standardErrors_mat = zeros(3,nSimu);
exitflag_mat = zeros(1,nSimu);

tic
for i = 1:nSimu
    [choices,iX] = simulateData(deltaU,capPi,nPeriods,nFirms);
    %{
    \subsection{Nested Fixed Point Maximum Likelihood Estimation}

    First, suppose that $\Pi$ is known. We use |fmincon| from \textsc{Matlab}'s \textsc{Optimization Toolbox} to maximize the partial likelihood for the choices (the code can easily be adapted to use other optimizers and packages, because these have a very similar \url{http://www.mathworks.nl/help/optim/ug/fmincon.html}{syntax}; see below). Because |fmincon| is a minimizer, we use minus the log likelihood as its objective. The function |negLogLik| computes this objective, but has input arguments other than the vector of model parameters to be estimated. Because \url{http://www.mathworks.nl/help/optim/ug/passing-extra-parameters.html}{the syntax of |fmincon| does not allow this}, we define a function handle |objectiveFunction| to an anonymous function that equals |negLogLik| but does not have this extra inputs.
    %}
    objectiveFunction = @(parameters)negLogLik(choices,iX,supportX,capPi,parameters(1:2),[delta(1);parameters(3)],...
                                               rho,@flowpayoffs,@bellman,@fixedPoint,tolFixedPoint)
    %{
    Before we can put |fmincon| to work on this objective function, we first have to set some of its other input arguments. We specify a $3\times 1$ vector |startvalues| with starting values for the parameters to be estimated, $(\beta_0,\beta_1,\delta_1)'$.
    %}
    startvalues = [-1;-0.1;0.5];
    %{
        We also set a lower bound of 0 on the third parameter, $\delta_1$, and (nonbinding) lower bounds of $-\infty$ on the other two parameters (|lowerBounds|). There is no need to specify upper bounds.\footnote{Note that |fmincon|, but also its alternatives discussed below, allow the user to specify bounds on parameters; if another function is used that does not allow for bounds on the parameters, you can use an alternative parameterization to ensure that parameters only take values in some admissible set (for example, you can specify $\delta_1=\exp(\delta_1^*)$ for $\delta_1^*\in\mathbb{R}$ to ensure that $\delta_1>0$). Minimizers like |fmincon| also allow you to impose more elaborate constraints on the parameters; you will need this option when implementing the MPEC alternative to NFXP of \cite{ecta12:juddsu} (see Section \ref{exercises}).}
    %}
    lowerBounds = -Inf*ones(size(startvalues));
    lowerBounds(3) = 0;
    %{
        Finally, we pass some options, including tolerances that specify the criterion for the outer loop convergence, to |fmincon| through the structure |OptimizerOptions| (recall that we have already set the inner loop tolerance in |tolFixedPoint|). We use the function |optimset| from the \textsc{Optimization Toolbox} to assign values to specific fields (options) in |OptimizerOptions| and then call |fmincon| to run the NFXP maximum likelihood procedure (to use \textsc{Knitro} instead, simply replace |fmincon| by |knitromatlab|, |knitrolink|, or |ktrlink|, depending on the packages installed\footnote{|fmincon| requires \textsc{Matlab}'s \textsc{Optimization Toolbox}, |knitromatlab| is included in \textsc{Knitro} 9.0, |knitrolink| uses both, and |ktrlink| can be used if the \textsc{Optimization Toolbox} is installed with an earlier version of \textsc{Knitro}.}).
    %}
    OptimizerOptions = optimset('Display','iter','Algorithm','interior-point','AlwaysHonorConstraints','bounds',...
                                'GradObj','on','TolFun',1E-6,'TolX',1E-10,'DerivativeCheck','off','TypicalX',[beta;delta(2)]);
    [maxLikEstimates,~,exitflag] = fmincon(objectiveFunction,startvalues,[],[],[],[],lowerBounds,[],[],OptimizerOptions)
    %{
    This gives maximum partial likelihood estimates of $(\beta_0,\beta_1,\delta_1)$. To calculate standard errors, we call |negLogLik| once more to estimate the corresponding Fisher information matrix and store this in |informationMatrix|. Its inverse is an estimate of the maximum likelihood estimator's asymptotic variance-covariance matrix.
    %}
    [~,~,informationMatrix] = objectiveFunction(maxLikEstimates);
    standardErrors = diag(sqrt(inv(informationMatrix)));
    %{
    The resulting parameter estimates and standard errors are displayed (third and fourth columns), together with the parameters' true (first column) and starting values (second column).
    %}
%     disp('Summary of Results');
%     disp('--------------------------------------------');
%     disp('      true     start     estim      ste.');
%     disp([[beta;delta(2)] startvalues maxLikEstimates standardErrors]);
    
    maxLikEstimates_mat(:,i) = maxLikEstimates;
    standardErrors_mat(:,i) = standardErrors;
    exitflag_mat(i) = exitflag;
end
toc

% disp(mean(maxLikEstimates_mat,2))
% disp(mean(standardErrors_mat,2))

disp('Summary of NFXP Monte Carlo Results');
disp('--------------------------------------------');
disp('      true     start   mean(estim)  ste.   mean(asymp ste.)');
disp([[beta;delta(2)] startvalues mean(maxLikEstimates_mat,2) std(maxLikEstimates_mat,0,2) mean(standardErrors_mat,2)]);

disp(['exitflags: ' num2str(exitflag_mat)])
disp(' ')
MSE = mean((maxLikEstimates_mat-repmat([beta;delta(2)],1,nSimu)).^2,2);
disp('MeanSquaredError')
disp(num2str(MSE));
%{
\subsection{Extension to an Unknown Markov Transition Matrix for the State}
Finally, consider the more realistic case that $\Pi$ is not known. In this case, \cite{nh94:rust} suggests a two-stage procedure. In the first stage, we estimate $\Pi$ using |estimatePi| and store the results in a $K\times K$ matrix |piHat|.
%}
piHat = estimatePi(iX,nSuppX)

NFXPmaxLikEstimates = maxLikEstimates_mat;
NFXPstandardErrors = standardErrors_mat;
NFXPexitflag = exitflag_mat;
NFXP_MSE = MSE;

%% MPEC
maxLikEstimates_mat = zeros(23,nSimu);
standardErrors_mat = zeros(23,nSimu);
exitflag_mat = zeros(1,nSimu);

tic
for i = 1:nSimu
    [choices,iX] = simulateData(deltaU,capPi,nPeriods,nFirms);
    %{
    \subsection{Nested Fixed Point Maximum Likelihood Estimation}

    First, suppose that $\Pi$ is known. We use |fmincon| from \textsc{Matlab}'s \textsc{Optimization Toolbox} to maximize the partial likelihood for the choices (the code can easily be adapted to use other optimizers and packages, because these have a very similar \url{http://www.mathworks.nl/help/optim/ug/fmincon.html}{syntax}; see below). Because |fmincon| is a minimizer, we use minus the log likelihood as its objective. The function |negLogLik| computes this objective, but has input arguments other than the vector of model parameters to be estimated. Because \url{http://www.mathworks.nl/help/optim/ug/passing-extra-parameters.html}{the syntax of |fmincon| does not allow this}, we define a function handle |objectiveFunction| to an anonymous function that equals |negLogLik| but does not have this extra inputs.
    %}
    % objectiveFunction = @(parameters)negLogLik(choices,iX,supportX,capPi,parameters(1:2),[delta(1);parameters(3)],...
    %                                            rho,@flowpayoffs,@bellman,@fixedPoint,tolFixedPoint)

    objectiveFunctionMpec = @(parameters)negLogLik_CE3(choices,iX,supportX,capPi,parameters(1:2),...
        [delta(1);parameters(3)],rho,@flowpayoffs,@bellman,@fixedPoint,tolFixedPoint,...
        parameters(4:13),parameters(14:23));

    %{
    Before we can put |fmincon| to work on this objective function, we first have to set some of its other input arguments. We specify a $3\times 1$ vector |startvalues| with starting values for the parameters to be estimated, $(\beta_0,\beta_1,\delta_1)'$.
    %}
    startvalues = [-1;-0.1;0.5;ones(20,1)];
    %{
        We also set a lower bound of 0 on the third parameter, $\delta_1$, and (nonbinding) lower bounds of $-\infty$ on the other two parameters (|lowerBounds|). There is no need to specify upper bounds.\footnote{Note that |fmincon|, but also its alternatives discussed below, allow the user to specify bounds on parameters; if another function is used that does not allow for bounds on the parameters, you can use an alternative parameterization to ensure that parameters only take values in some admissible set (for example, you can specify $\delta_1=\exp(\delta_1^*)$ for $\delta_1^*\in\mathbb{R}$ to ensure that $\delta_1>0$). Minimizers like |fmincon| also allow you to impose more elaborate constraints on the parameters; you will need this option when implementing the MPEC alternative to NFXP of \cite{ecta12:juddsu} (see Section \ref{exercises}).}
    %}
    lowerBounds = -Inf*ones(size(startvalues));
    lowerBounds(3) = 0;
    %{
        Finally, we pass some options, including tolerances that specify the criterion for the outer loop convergence, to |fmincon| through the structure |OptimizerOptions| (recall that we have already set the inner loop tolerance in |tolFixedPoint|). We use the function |optimset| from the \textsc{Optimization Toolbox} to assign values to specific fields (options) in |OptimizerOptions| and then call |fmincon| to run the NFXP maximum likelihood procedure (to use \textsc{Knitro} instead, simply replace |fmincon| by |knitromatlab|, |knitrolink|, or |ktrlink|, depending on the packages installed\footnote{|fmincon| requires \textsc{Matlab}'s \textsc{Optimization Toolbox}, |knitromatlab| is included in \textsc{Knitro} 9.0, |knitrolink| uses both, and |ktrlink| can be used if the \textsc{Optimization Toolbox} is installed with an earlier version of \textsc{Knitro}.}).
    %}
    % OptimizerOptions = optimset('Display','iter','Algorithm','interior-point','AlwaysHonorConstraints','bounds',...
    %                             'GradObj','on','TolFun',1E-6,'TolX',1E-10,'DerivativeCheck','off','TypicalX',[beta;delta(2)]);
    OptimizerOptions = optimset('Display','iter','Algorithm','interior-point','AlwaysHonorConstraints','bounds',...
                            'GradObj','off','TolFun',1E-6,'TolX',1E-10,'DerivativeCheck','off', 'MaxFunEvals',30000);
    [maxLikEstimates,~,exitflag,~,~,~,hessian] = fmincon(objectiveFunctionMpec,startvalues,[],[],[],[],lowerBounds,[],@(parameters)(constraintMpec(parameters,delta, rho, capPi, supportX)),OptimizerOptions)
    %{
    This gives maximum partial likelihood estimates of $(\beta_0,\beta_1,\delta_1)$. To calculate standard errors, we call |negLogLik| once more to estimate the corresponding Fisher information matrix and store this in |informationMatrix|. Its inverse is an estimate of the maximum likelihood estimator's asymptotic variance-covariance matrix.
    %}
    % [~,~,informationMatrix] = objectiveFunction(maxLikEstimates);        % tbf
    [~,~,informationMatrix] = objectiveFunctionMpec(maxLikEstimates);        % tbf
    % standardErrors = diag(sqrt(inv(informationMatrix)));
    standardErrors = sqrt(diag(hessian).^-1);
    
    maxLikEstimates_mat(:,i) = maxLikEstimates;
    standardErrors_mat(:,i) = standardErrors;
    exitflag_mat(i) = exitflag;
end
toc

% disp(mean(maxLikEstimates_mat,2))
% disp(mean(standardErrors_mat,2))

disp('Summary of MPEC Monte Carlo Results');
disp('--------------------------------------------');
disp('      true     start   mean(estim)  ste.   mean(asymp ste.)');
disp([[beta;delta(2)] startvalues(1:3) mean(maxLikEstimates_mat(1:3,:),2) std(maxLikEstimates_mat(1:3,:),0,2) mean(standardErrors_mat(1:3,:),2)]);

% disp([[beta;delta(2)] startvalues(1:3) maxLikEstimates(1:3) standardErrors(1:3)]);

disp(['exitflags: ' num2str(exitflag_mat)])
disp(' ')
MSE = mean((maxLikEstimates_mat(1:3,:)-repmat([beta;delta(2)],1,nSimu)).^2,2);
disp('MeanSquaredError')
disp(num2str(MSE));

MPECmaxLikEstimates = maxLikEstimates_mat;
MPECstandardErrors = standardErrors_mat;
MPECexitflag = exitflag_mat;
MPEC_MSE = MSE;

%% CCP(HMSS)
maxLikEstimates_mat = zeros(3,nSimu);
standardErrors_mat = zeros(3,nSimu);
exitflag_mat = zeros(1,nSimu);

tic
for i = 1:nSimu
    [choices,iX] = simulateData(deltaU,capPi,nPeriods,nFirms);

    % compute estimated Delta U from nonparametetrically estimated chioce probabilities
    DeltaU = estimateDeltaU(choices,iX,nSuppX,nFirms,nPeriods);

    % draw Xs and epsilons to be used for forward simultation
    nSimuPeriods = 200;               % set how many periods to forward simulate
    nSimuFirms = 1000;                % set how many times to simulate for each starting point 
    oneMinPi = eye(nSuppX)-capPi';
    pInf     = [oneMinPi(1:nSuppX-1,:);ones(1,nSuppX)]\[zeros(nSuppX-1,1);1];
    drawsX = zeros(nSimuPeriods,nSimuFirms,nSuppX);
    % for k = 1:(nSuppX*2)
    %     drawsX = randomDiscrete(pInf*ones(1,nFirms));
    %     drawsX = zeros(1,nFirms);                     % Xt-1
    for k = 1:nSuppX
        drawsX_k = k.*ones(1,nSimuFirms);               % row 1 is Xt       Xt=1 here
        drawsEpsilon1 = random('ev',zeros(1,nSimuFirms),ones(1,nSimuFirms));
        drawsEpsilon0 = random('ev',zeros(1,nSimuFirms),ones(1,nSimuFirms));
        for t = 2:nSimuPeriods
            drawsX_k = [drawsX_k;randomDiscrete(capPi(drawsX_k(end,:),:)')];           % Xtau
            drawsEpsilon1 = [drawsEpsilon1;random('ev',zeros(1,nSimuFirms),ones(1,nSimuFirms))];
            drawsEpsilon0 = [drawsEpsilon0;random('ev',zeros(1,nSimuFirms),ones(1,nSimuFirms))];
        end
        drawsX(:,:,k) = drawsX_k;
    end

    %{
    \subsection{Nested Fixed Point Maximum Likelihood Estimation}

    First, suppose that $\Pi$ is known. We use |fmincon| from \textsc{Matlab}'s \textsc{Optimization Toolbox} to maximize the partial likelihood for the choices (the code can easily be adapted to use other optimizers and packages, because these have a very similar \url{http://www.mathworks.nl/help/optim/ug/fmincon.html}{syntax}; see below). Because |fmincon| is a minimizer, we use minus the log likelihood as its objective. The function |negLogLik| computes this objective, but has input arguments other than the vector of model parameters to be estimated. Because \url{http://www.mathworks.nl/help/optim/ug/passing-extra-parameters.html}{the syntax of |fmincon| does not allow this}, we define a function handle |objectiveFunction| to an anonymous function that equals |negLogLik| but does not have this extra inputs.
    %}
    objectiveFunction = @(parameters)negLogLik_CE5(choices,iX,supportX,capPi,parameters(1:2),[delta(1);parameters(3)],...
                                               rho,drawsX,drawsEpsilon1,drawsEpsilon0,DeltaU,...
                                               @flowpayoffs,@bellman,@fixedPoint,tolFixedPoint)
    %{
    Before we can put |fmincon| to work on this objective function, we first have to set some of its other input arguments. We specify a $3\times 1$ vector |startvalues| with starting values for the parameters to be estimated, $(\beta_0,\beta_1,\delta_1)'$.
    %}
    startvalues = [-1;-0.1;0.5];
    %{
        We also set a lower bound of 0 on the third parameter, $\delta_1$, and (nonbinding) lower bounds of $-\infty$ on the other two parameters (|lowerBounds|). There is no need to specify upper bounds.\footnote{Note that |fmincon|, but also its alternatives discussed below, allow the user to specify bounds on parameters; if another function is used that does not allow for bounds on the parameters, you can use an alternative parameterization to ensure that parameters only take values in some admissible set (for example, you can specify $\delta_1=\exp(\delta_1^*)$ for $\delta_1^*\in\mathbb{R}$ to ensure that $\delta_1>0$). Minimizers like |fmincon| also allow you to impose more elaborate constraints on the parameters; you will need this option when implementing the MPEC alternative to NFXP of \cite{ecta12:juddsu} (see Section \ref{exercises}).}
    %}
    lowerBounds = -Inf*ones(size(startvalues));
    lowerBounds(3) = 0;
    %{
        Finally, we pass some options, including tolerances that specify the criterion for the outer loop convergence, to |fmincon| through the structure |OptimizerOptions| (recall that we have already set the inner loop tolerance in |tolFixedPoint|). We use the function |optimset| from the \textsc{Optimization Toolbox} to assign values to specific fields (options) in |OptimizerOptions| and then call |fmincon| to run the NFXP maximum likelihood procedure (to use \textsc{Knitro} instead, simply replace |fmincon| by |knitromatlab|, |knitrolink|, or |ktrlink|, depending on the packages installed\footnote{|fmincon| requires \textsc{Matlab}'s \textsc{Optimization Toolbox}, |knitromatlab| is included in \textsc{Knitro} 9.0, |knitrolink| uses both, and |ktrlink| can be used if the \textsc{Optimization Toolbox} is installed with an earlier version of \textsc{Knitro}.}).
    %}
    OptimizerOptions = optimset('Display','iter','Algorithm','interior-point','AlwaysHonorConstraints','bounds',...
                                'GradObj','off','TolFun',1E-6,'TolX',1E-10,'DerivativeCheck','off');
    [maxLikEstimates,~,exitflag,~,~,~,hessian] = fmincon(objectiveFunction,startvalues,[],[],[],[],lowerBounds,[],[],OptimizerOptions)
    %{
    This gives maximum partial likelihood estimates of $(\beta_0,\beta_1,\delta_1)$. To calculate standard errors, we call |negLogLik| once more to estimate the corresponding Fisher information matrix and store this in |informationMatrix|. Its inverse is an estimate of the maximum likelihood estimator's asymptotic variance-covariance matrix.
    %}
    % [~,~,informationMatrix] = objectiveFunction(maxLikEstimates);

    standardErrors = diag(sqrt(inv(hessian)));
    
    maxLikEstimates_mat(:,i) = maxLikEstimates;
    standardErrors_mat(:,i) = standardErrors;
    exitflag_mat(i) = exitflag;
end
toc

% disp(mean(maxLikEstimates_mat,2))
% disp(mean(standardErrors_mat,2))

disp('Summary of CCP Monte Carlo Results');
disp('--------------------------------------------');
disp('      true     start   mean(estim)  ste.   mean(asymp ste.)');
disp([[beta;delta(2)] startvalues mean(maxLikEstimates_mat,2) std(maxLikEstimates_mat,0,2) mean(standardErrors_mat,2)]);

disp(['exitflags: ' num2str(exitflag_mat)])
disp(' ')
MSE = mean((maxLikEstimates_mat-repmat([beta;delta(2)],1,nSimu)).^2,2);
disp('MeanSquaredError')
disp(num2str(MSE));

CCPmaxLikEstimates = maxLikEstimates_mat;
CCPstandardErrors = standardErrors_mat;
CCPexitflag = exitflag_mat;
CCP_MSE = MSE;
