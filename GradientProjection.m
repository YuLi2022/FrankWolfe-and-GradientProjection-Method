function [xk, x_iters]= GradientProjection(x0, func, grad_func, A, b, Aeq, beq,...
    eps, maxIter, showInfo) 
% 梯度投影法
% 函数说明：均要写成这个形式
%   min func(X)
%   s.t. A*X >= b; Aeq*X == beq; 
% 输入:
%   x0: 初始迭代解；
%   func, grad_func: 目标函数及其梯度函数的句柄
%   eps: 求解精度
%   maxIter: 允许的最大迭代次数，不能无限循环 
%   showInfo: true输出迭代信息；false不输出迭代信息 
% 输出：
%   xk: 最后一次迭代解, 相当于是找到的最优解  
%   x_iters: 保存每一次迭代的解 
findKT = false; 

xk = x0; 
k = 1; 
x_iters = xk'; 
if showInfo
    disp("当前正在第" + num2str(k) + "次迭代...");
    disp("x(" + num2str(k) + ") = "); 
    disp(num2str(xk)); 
end
% 首先在xk处分解为A1,b1和A2,b2 
[A1, b1, A2, b2] = split(A, b, xk, 1e-4); % 步骤2 

if isempty(maxIter)  
    maxIter = 1e3; 
end
while (~findKT) && (k < maxIter)
    % 步骤3，构造M，计算P
    M = [A1;Aeq]; 
    if isempty(M) 
        % 如果M为空
        P = eye(size(A,2));  % 构造单位矩阵 
    else
        P = eye(size(A,2)) - M' * inv(M*M') * M; 
    end
    
    % 步骤4
    dk = -P*grad_func(xk);  
    if all(abs(dk)<=eps) 
        % 如果d(k)=0, 
        % 步骤5
        if isempty(M) 
            findKT=true;  % 算法停止计算 
        else
            W = inv(M*M')*M*grad_func(xk);
            nU = size(A1,1);  % 
            u = W(1:nU); 
            if all(u>=0)  
                % 停止计算，当前解为KT点
                findKT=true; 
            else
                % 选择一个负分量，修正A1，去掉A1中对应的uj的行，重新计算P
                indices = 1:nU; 
                ids = indices(u<0);  % 小于0的索引
                [~, idx] = min(u(u<0));  % 选择相应的索引 
                idx = ids(idx); 
                A1(idx,:) = [];  % 去掉这一行
                continue; 
            end
        end
    else
        % 步骤6； 
        b_ = b2 - A2*xk; 
        d_ = A2*dk; 
        if all(d_>=0)  
            lambda_max = 10; 
        else
            lambda_max = min(b_(d_<0)./d_(d_<0)); 
        end
        lambda = golden_selection(0, lambda_max, 1e-6, xk, dk, func);
        xk = xk + lambda*dk; 
        k = k + 1;  % 迭代次数+1 
        [A1, b1, A2, b2] = split(A, b, xk, 1e-4); % 步骤2 
        x_iters = [x_iters;xk']; 
        if showInfo
            disp("当前正在第" + num2str(k) + "次迭代...");
            disp("x(" + num2str(k) + ") = "); 
            disp(num2str(xk))
        end
    end
end
end


function [A1, b1, A2, b2] = split(A, b, xk, eps)   
% eps: 判断精度 
isch = abs(A*xk - b) <= eps;  
indices = 1:size(A,1); 
A1 = A(indices(isch),:); 
b1 = b(indices(isch)); 
A2 = A(indices(~isch),:); 
b2 = b(indices(~isch)); 
end 

function ak = golden_selection(a, b, eps, xk, dk, func) 
% 求解当前最优步长
% [a,b]: 搜索区间 
% eps: 精确度 
t = 0.618;  
ak = a; 
bk = b; 
while (bk-ak) > eps
    lk = ak + (1-t)*(bk-ak); 
    rk = ak + t*(bk-ak); 
    if func(xk+lk*dk) > func(xk+rk*dk) 
        ak = lk; 
    else
        bk = rk; 
    end
end
end