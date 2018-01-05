function [LL]=SML_ObjectFun(b0)               %%% w: also need to restrict SEs to be above zeros ! tbc

global DATA1 DATA2 v_IF rep_v_IF NF NI_F NS NJ NObs NI NT 

alpha = b0(1);
beta = b0(2);
mu = b0(3);
sigma = b0(4);
delta = b0(5);

%% compute the demand probabilities
u_ifjt = repmat(DATA1(:,5:6)*[beta alpha]', [1, NS])+mu+ sigma*rep_v_IF;       % NObs*NS      75000*200
%%% NP: This repmat can be done in the Main file, so you don't need to
%%% repeat it at each iteration of the objective function

u_ifjt(DATA1(:,4)==6,:)=0;                     % set the outside option utility to be 0

% compute the sum of exp(utility) for all 6 options
year_index = (1:NObs/(NJ+1))*6;
exp_u_sum_j = zeros(NObs,NS);                                                   % NObs*NS      75000*200
for i = 1:length(year_index)
    u_year_t = u_ifjt((year_index(i)-5):year_index(i),:);                       % (NJ+1)*NS       5*200
    temp_sum = sum(exp(u_year_t),1);                                            % 1*NS        1*200
    exp_u_sum_j((year_index(i)-5):year_index(i),:) = repmat(temp_sum,NJ+1,1);   % (NJ+1)*NS       5*200
end
%%% NP: This is a short loop so it's fine, but it could be vectorized to
%%% improve speed
 
% compute probability of choosing each option
 PD_ifjt = exp(u_ifjt)./exp_u_sum_j;                         % NObs*NS      75000*200
 chosen_PD_ifjt = PD_ifjt.*repmat(DATA1(:,7),1,NS);          % NObs*NS      75000*200
 chosen_PD_ifjt(chosen_PD_ifjt==0)=1;                        % avoid numerical issue with log by setting 0 values to 1
 log_chosen_PD_ifjt = log(chosen_PD_ifjt);
 
 % summing up the log probabilities 
 LL_D = 0;
 row_index = (1:NObs/30)*30;                        % 1*2500
 for i = 1:length(row_index)                        % summing over log probability of each individual i working in f firms
     sum_log_chosen_PD_ifjt = sum(log_chosen_PD_ifjt(row_index(i)-29:row_index(i),:),1);   % 1*NS  1*200         % only summing the *chosen* log probabilities of each individual i on year t
     %%% NP: A faster way of considering only the "chosen" probabilities is
     %%% calculating the product of (Prob.^Data)
     temp_PD = exp(sum_log_chosen_PD_ifjt);
     sum_s_PD = mean(temp_PD,2);                    % 1*1                    % average over NS simulation draws
     log_sum_s_PD = log(sum_s_PD);
     LL_D = LL_D + log_sum_s_PD;                   
 end
 %%% NP: This loop could be avoided by arranging your data in Main such
 %%% that the probabilities are in a 5-dimensional matrix of dimensions
 %%% (years,individual,product,firm,simulation), then you can calculate
 %%% the demand log-likelihood as:
 %%%  llike = -sum(sum(log(mean(prod(prod(Prob.^Choice,3),1),5)),2),4);
 
 
 %% compute the cost probabilities
 % computing expected risk of individuals working in each firm f
 fj_index = (1:5)*500;                            % 1*5
 rep_index = (1:5)*25;                            % 1*5
 rep_risk_fj = zeros(size(DATA2,1),NS);           % 125*200
 for i= 1:length(fj_index)
     mean_risk_f = mean(v_IF(fj_index(i)-499:fj_index(i),:),1);    % 1*200
     temp_risk_fj = repmat(mean_risk_f,25,1);                      % 25*200
     rep_risk_fj(rep_index(i)-24:rep_index(i),:) = temp_risk_fj;   % replicate risks to have the same # of rows as in DATA2
 end
 %%% NP: This loop could be also avoided following the previous comment
 
% compute probability of observing the cost residuals
rep_C_fjt = repmat(DATA2(:,4),1,NS);             % 125*200
PC_fjts = normpdf(rep_C_fjt-delta-sigma.*rep_risk_fj);         % 125*200
 
% test
% A = magic(6)
% AA =reshape(A',[6,2,3])
% permute(AA,[2 1 3])

% preparing the probability matrix for taking average over simulation draws
PC_fjts_temp = reshape(PC_fjts',[NS,5,25]);
PC_fjts_cube = permute(PC_fjts_temp,[2 1 3]);     % NJ*NS*(NFxNT) 5*200*25   
log_PC_fjts_cube = log(PC_fjts_cube);

t_index = (1:5)*5;                             % 1*5
PC_fjs_cube = zeros(5,NS,5);                   % NJ*NS*NF        5*200*5
for i = 1:length(t_index)
    PC_fjs_cube(:,:,i)= exp(sum(log_PC_fjts_cube(:,:,t_index(i)-4:t_index(i)),3));
end

% taking average over simulation draws
PC_fjs_cube_temp = permute(PC_fjs_cube,[1 3 2]);    % NJ*NF*NS        5*5*200
PC_fj = mean(PC_fjs_cube_temp,3);                   % NJ*NF           5*5           

% summing up log cost probabilities
LL_C = sum(sum(log(PC_fj)));                        % 1*1

% combining demand and cost probabilities
LL = -(LL_D + LL_C);
    
    

 

    
