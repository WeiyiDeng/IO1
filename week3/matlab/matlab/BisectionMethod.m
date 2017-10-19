%% bisection method_numerical method for solving nonlinear equation
% w Oct 17 2017
tic
tryfun = @(x) 0.6*x+1000*log(1+0.02*x)-3000;
x_L = 0;
x_U = 1000;
x_mid = 500;
while abs(tryfun(x_mid)) > eps
    if tryfun(x_L)*tryfun(x_mid)<0
        x_U = x_mid;
    elseif tryfun(x_U)*tryfun(x_mid)<0
        x_L = x_mid;
    end
    x_mid = (x_L+x_U)/2;
end
x_mid
toc

% compare to Zeroin Algorithm by matlab
x = fzero(tryfun,[0 1000],optimset('display','iter'));     
tic
x = fzero(tryfun,[0 1000]);
x
toc