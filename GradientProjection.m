function [xk, x_iters]= GradientProjection(x0, func, grad_func, A, b, Aeq, beq,...
    eps, maxIter, showInfo) 
% �ݶ�ͶӰ��
% ����˵������Ҫд�������ʽ
%   min func(X)
%   s.t. A*X >= b; Aeq*X == beq; 
% ����:
%   x0: ��ʼ�����⣻
%   func, grad_func: Ŀ�꺯�������ݶȺ����ľ��
%   eps: ��⾫��
%   maxIter: �������������������������ѭ�� 
%   showInfo: true���������Ϣ��false�����������Ϣ 
% �����
%   xk: ���һ�ε�����, �൱�����ҵ������Ž�  
%   x_iters: ����ÿһ�ε����Ľ� 
findKT = false; 

xk = x0; 
k = 1; 
x_iters = xk'; 
if showInfo
    disp("��ǰ���ڵ�" + num2str(k) + "�ε���...");
    disp("x(" + num2str(k) + ") = "); 
    disp(num2str(xk)); 
end
% ������xk���ֽ�ΪA1,b1��A2,b2 
[A1, b1, A2, b2] = split(A, b, xk, 1e-4); % ����2 

if isempty(maxIter)  
    maxIter = 1e3; 
end
while (~findKT) && (k < maxIter)
    % ����3������M������P
    M = [A1;Aeq]; 
    if isempty(M) 
        % ���MΪ��
        P = eye(size(A,2));  % ���쵥λ���� 
    else
        P = eye(size(A,2)) - M' * inv(M*M') * M; 
    end
    
    % ����4
    dk = -P*grad_func(xk);  
    if all(abs(dk)<=eps) 
        % ���d(k)=0, 
        % ����5
        if isempty(M) 
            findKT=true;  % �㷨ֹͣ���� 
        else
            W = inv(M*M')*M*grad_func(xk);
            nU = size(A1,1);  % 
            u = W(1:nU); 
            if all(u>=0)  
                % ֹͣ���㣬��ǰ��ΪKT��
                findKT=true; 
            else
                % ѡ��һ��������������A1��ȥ��A1�ж�Ӧ��uj���У����¼���P
                indices = 1:nU; 
                ids = indices(u<0);  % С��0������
                [~, idx] = min(u(u<0));  % ѡ����Ӧ������ 
                idx = ids(idx); 
                A1(idx,:) = [];  % ȥ����һ��
                continue; 
            end
        end
    else
        % ����6�� 
        b_ = b2 - A2*xk; 
        d_ = A2*dk; 
        if all(d_>=0)  
            lambda_max = 10; 
        else
            lambda_max = min(b_(d_<0)./d_(d_<0)); 
        end
        lambda = golden_selection(0, lambda_max, 1e-6, xk, dk, func);
        xk = xk + lambda*dk; 
        k = k + 1;  % ��������+1 
        [A1, b1, A2, b2] = split(A, b, xk, 1e-4); % ����2 
        x_iters = [x_iters;xk']; 
        if showInfo
            disp("��ǰ���ڵ�" + num2str(k) + "�ε���...");
            disp("x(" + num2str(k) + ") = "); 
            disp(num2str(xk))
        end
    end
end
end


function [A1, b1, A2, b2] = split(A, b, xk, eps)   
% eps: �жϾ��� 
isch = abs(A*xk - b) <= eps;  
indices = 1:size(A,1); 
A1 = A(indices(isch),:); 
b1 = b(indices(isch)); 
A2 = A(indices(~isch),:); 
b2 = b(indices(~isch)); 
end 

function ak = golden_selection(a, b, eps, xk, dk, func) 
% ��⵱ǰ���Ų���
% [a,b]: �������� 
% eps: ��ȷ�� 
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