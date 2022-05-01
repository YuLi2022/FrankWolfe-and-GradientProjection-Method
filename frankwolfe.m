function  [xk, x_iters] = frankwolfe(x0, func, grad_func, A, b, Aeq, beq, lb, ub, ...
    eps, maxIter, showInfo) 
% frank worf算法
% 输入优化问题示例：
%   min func
%   s.t. A*X<=b; Aeq*X=beq;  lb<= X <= ub; 
% 输入:
%   x0: 初始迭代解；
%   func, grad_func: 目标函数及其梯度函数的句柄
%   eps: 求解精度
%   maxIter: 允许的最大迭代次数，不能无限循环 
%   showInfo: true输出迭代信息；false不输出迭代信息 
% 输出：
%   xk: 最后一次迭代解, 相当于是找到的最优解  
%   x_iters: 保存每一次迭代的解 

xk = x0; 
grad_xk = grad_func(xk);  % 计算梯度; 
options = optimoptions('linprog', 'Display','off'); 
yk = linprog(grad_xk,A,b,Aeq,beq,lb,ub,options);

% 记录信息 
x_iters = xk'; 
k = 1; 
if showInfo
    disp("当前正在第" + num2str(k) + "次迭代...");
    disp("x(" + num2str(k) + ") = " + num2str(xk)); 
end
if isempty(maxIter)
    maxIter = 1e3;  % 最大搜索次数 
end

while (abs(grad_xk'*(yk-xk)) >= eps) && (k<maxIter) 
    % 进行一维搜索，
    lambda = golden_selection(0, 1, 1e-6, xk, yk, func); 
    xk = xk + lambda*(yk-xk);   % 新的x
    grad_xk = grad_func(xk);  % 计算梯度; 
    yk = linprog(grad_xk,A,b,Aeq,beq,lb,ub,options); 
    % 保存信息 
    x_iters = [x_iters;xk']; 
    k = k + 1;
    if showInfo
        disp("当前正在第" + num2str(k) + "次迭代...");
        disp("x(" + num2str(k) + ") = " + num2str(xk));
    end
end
end

function ak = golden_selection(a, b, eps, xk, yk, func) 
% 求解当前最优步长
% [a,b]: 搜索区间 
% eps: 精确度 
t = 0.618;  
ak = a; 
bk = b; 
while (bk-ak) > eps
    lk = ak + (1-t)*(bk-ak); 
    rk = ak + t*(bk-ak); 
    if func(xk+lk*(yk-xk)) > func(xk+rk*(yk-xk)) 
        ak = lk; 
    else
        bk = rk; 
    end
end
end