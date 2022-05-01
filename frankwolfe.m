function  [xk, x_iters] = frankwolfe(x0, func, grad_func, A, b, Aeq, beq, lb, ub, ...
    eps, maxIter, showInfo) 
% frank worf�㷨
% �����Ż�����ʾ����
%   min func
%   s.t. A*X<=b; Aeq*X=beq;  lb<= X <= ub; 
% ����:
%   x0: ��ʼ�����⣻
%   func, grad_func: Ŀ�꺯�������ݶȺ����ľ��
%   eps: ��⾫��
%   maxIter: �������������������������ѭ�� 
%   showInfo: true���������Ϣ��false�����������Ϣ 
% �����
%   xk: ���һ�ε�����, �൱�����ҵ������Ž�  
%   x_iters: ����ÿһ�ε����Ľ� 

xk = x0; 
grad_xk = grad_func(xk);  % �����ݶ�; 
options = optimoptions('linprog', 'Display','off'); 
yk = linprog(grad_xk,A,b,Aeq,beq,lb,ub,options);

% ��¼��Ϣ 
x_iters = xk'; 
k = 1; 
if showInfo
    disp("��ǰ���ڵ�" + num2str(k) + "�ε���...");
    disp("x(" + num2str(k) + ") = " + num2str(xk)); 
end
if isempty(maxIter)
    maxIter = 1e3;  % ����������� 
end

while (abs(grad_xk'*(yk-xk)) >= eps) && (k<maxIter) 
    % ����һά������
    lambda = golden_selection(0, 1, 1e-6, xk, yk, func); 
    xk = xk + lambda*(yk-xk);   % �µ�x
    grad_xk = grad_func(xk);  % �����ݶ�; 
    yk = linprog(grad_xk,A,b,Aeq,beq,lb,ub,options); 
    % ������Ϣ 
    x_iters = [x_iters;xk']; 
    k = k + 1;
    if showInfo
        disp("��ǰ���ڵ�" + num2str(k) + "�ε���...");
        disp("x(" + num2str(k) + ") = " + num2str(xk));
    end
end
end

function ak = golden_selection(a, b, eps, xk, yk, func) 
% ��⵱ǰ���Ų���
% [a,b]: �������� 
% eps: ��ȷ�� 
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