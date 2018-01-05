%%
% Xt = 1;
% At_prev = 0;
% A_hat_t = DeltaU(Xt,At_prev+1)+drawsEpsilon1(1,:)-drawsEpsilon0(1,:)>=0;
% A_hat = zeros(size(drawsX));
% for t = 2:nSimuPeriods
%     X_hat_tau = drawsX(t,:);            % tau = t+1
%     Epsilon0_hat_tau = drawsEpsilon0(t,:);
%     Epsilon1_hat_tau = drawsEpsilon1(t,:);
%     A_hat_tau = diag(DeltaU(X_hat_tau,A_hat_t+1))' + Epsilon1_hat_tau - Epsilon0_hat_tau>=0;
%     A_hat(t,:) = A_hat_tau;
%     A_hat_t = A_hat_tau;
% end

% debug
nFirms = 1000;
nSimuPeriods = 1000;                                  %
oneMinPi = eye(nSuppX)-capPi';
pInf     = [oneMinPi(1:nSuppX-1,:);ones(1,nSuppX)]\[zeros(nSuppX-1,1);1];
drawsX = zeros(nSimuPeriods,nFirms,nSuppX);
drawsEpsilon1 = zeros(nSimuPeriods,nFirms,nSuppX);
drawsEpsilon0 = zeros(nSimuPeriods,nFirms,nSuppX);
% for k = 1:(nSuppX*2)
%     drawsX = randomDiscrete(pInf*ones(1,nFirms));
%     drawsX = zeros(1,nFirms);                     % Xt-1
for k = 1:nSuppX
    drawsX_k = k.*ones(1,nFirms);            % row 1 is Xt       Xt=1 here
    drawsEpsilon1_k = random('ev',zeros(1,nFirms),ones(1,nFirms));
    drawsEpsilon0_k = random('ev',zeros(1,nFirms),ones(1,nFirms));
%     drawsX_k = repmat(drawsX_k,nSimuPeriods,1);    
%     drawsEpsilon1_k = repmat(drawsEpsilon1_k,nSimuPeriods,1);
%     drawsEpsilon0_k = repmat(drawsEpsilon0_k,nSimuPeriods,1);
    for t = 2:nSimuPeriods
        drawsX_k = [drawsX_k;randomDiscrete(capPi(drawsX_k(end,:),:)')];           % Xtau
        drawsEpsilon1_k = [drawsEpsilon1_k;random('ev',zeros(1,nFirms),ones(1,nFirms))];
        drawsEpsilon0_k = [drawsEpsilon0_k;random('ev',zeros(1,nFirms),ones(1,nFirms))];
    end
    drawsX(:,:,k) = drawsX_k;
    drawsEpsilon1(:,:,k) = drawsEpsilon1_k;
    drawsEpsilon0(:,:,k) = drawsEpsilon0_k;
end
DeltaU = deltaU;                                % (Xt, At-1)

% euler_constant = double(eulergamma);
euler_constant = 0;
%
[u0,u1] = flowpayoffs(supportX,beta,delta);         
u_cube = u0;
u_cube(:,:,2) = u1;
nSuppX = size(supportX,1);
A_hat_cell = cell(nSuppX,2);
u_hat_cell = cell(nSuppX,2);
for Xt = 1:nSuppX
% for Xt = 1
    for At_prev = [0 1]
%     for At_prev = 0                             % At-1
        A_hat_t = DeltaU(Xt,At_prev+1)+drawsEpsilon1(1,:,Xt)-drawsEpsilon0(1,:,Xt)-euler_constant>=0;
        u_hat_At = zeros(size(A_hat_t));
        for i = 1:length(A_hat_t)
            u_hat_At(i) = u_cube(Xt,At_prev+1,A_hat_t(i)+1);
        end
        A_hat = zeros(nSimuPeriods,nFirms);
        u_hat = zeros(nSimuPeriods,nFirms);
        A_hat(1,:) = A_hat_t;
        u_hat(1,:) = u_hat_At;
        for t = 2:nSimuPeriods
            X_hat_tau = drawsX(t,:,Xt);            % tau = t+1
            Epsilon0_hat_tau = drawsEpsilon0(t,:,Xt);
            Epsilon1_hat_tau = drawsEpsilon1(t,:,Xt);
            A_hat_tau = zeros(size(A_hat_t));
            for i = 1:length(A_hat_t)
                A_hat_tau(i) = DeltaU(X_hat_tau(i),A_hat_t(i)+1) + Epsilon1_hat_tau(i) - Epsilon0_hat_tau(i)-euler_constant>=0;
            end
%             A_hat_tau = diag(DeltaU(X_hat_tau,A_hat_t+1))' + Epsilon1_hat_tau - Epsilon0_hat_tau>=0;
            A_hat(t,:) = A_hat_tau;
            for i = 1:length(A_hat_t)
                u_hat(t,i) = u_cube(X_hat_tau(i),A_hat_t(i)+1,A_hat_tau(i)+1);
            end
            A_hat_t = A_hat_tau;
        end
        A_hat_cell{Xt,At_prev+1} = A_hat;
        u_hat_cell{Xt,At_prev+1} = u_hat;
    end
end

CapU_cube_hat = zeros(nSuppX,2,2);
Epsilon_quad = drawsEpsilon0;
Epsilon_quad(:,:,:,2) = drawsEpsilon1;
for Xt = 1:nSuppX
    for At_prev = [0 1]
        A_hat11 = A_hat_cell{Xt,At_prev+1};
        u_hat11 = u_hat_cell{Xt,At_prev+1};
        CapU11 = u_hat11(1,:);
        Epsilon_Atau = zeros(size(A_hat11));
        for t = 2:nSimuPeriods
            for i = 1:nFirms
                Epsilon_Atau(t,i) = Epsilon_quad(t,i,Xt,A_hat11(t,i)+1);
            end
            CapU11 = CapU11+rho^(t-1)*(u_hat11(t,:)+Epsilon_Atau(t,:));
        end
        CapU11_bar = mean(CapU11);
        CapU11_bar_At0 = mean(CapU11(A_hat11(1,:)==0));
        CapU11_bar_At1 = mean(CapU11(A_hat11(1,:)==1));
        delta_CapU11 = CapU11_bar_At1-CapU11_bar_At0;
        CapU_cube_hat(Xt,At_prev+1,1) = CapU11_bar_At0;
        CapU_cube_hat(Xt,At_prev+1,2) = CapU11_bar_At1;
    end
end
delta_CapU_hat = CapU_cube_hat(:,:,2)- CapU_cube_hat(:,:,1)