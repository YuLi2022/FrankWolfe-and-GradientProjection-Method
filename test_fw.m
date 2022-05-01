% % frank worf�㷨�������� 
% clear; 
% func = @objective; 
% grad_func = @grad_objective; 
% 
% A = []; 
% b = []; 
% Aeq = []; 
% beq = []; 
% lb = [0;0]; 
% ub = [5;8]; 
% x0=[0;0]; 
% eps = 1e-3; 
% maxIter = 2e3; % ����������100�Σ� 
% % ����frankworf�㷨 
% [sol, x_iters] = frankwolfe(x0,func,grad_func,A,b,Aeq,beq,lb,ub,eps,maxIter,false); 
% disp("��������: "+ num2str(size(x_iters,1))); 
% disp("�ҵ������Ž�Ϊ:"); 
% disp(num2str(sol)); 
% 
% % ����ÿһ�ε�����ý�ĺ���ֵ 
% y_iters = zeros(1,size(x_iters,1)); 
% for i = 1:length(y_iters) 
%     y_iters(i) = func(x_iters(i,:)); 
% end
% figure();
% plot(1:length(y_iters), y_iters); 
% xlabel('��������');
% ylabel('Ŀ�꺯��ֵ'); 
% 
% % ����Ŀ�꺯�����ݶȺ��� 
% function y = objective(x) 
% y = 2*x(1)^2 + 3*x(2)^2 - 4*x(1)*x(2)-10*x(1) + 7; 
% end
% 
% function gx = grad_objective(x) 
% gx = [4*x(1)-4*x(2)-10; 6*x(2)-4*x(1)]; 
% end



func = @objective; 
grad_func = @grad_objective; 

% ���Ժ���  
A = []; 
b = []; 
Aeq = [1,1,1,0;1,5,0,1]; 
beq = [3;6]; 
lb = [0,0,0,0];  
ub = []; 
%x0 = [2;0;1;4];
x0 = [0.5;1;0.5;0.5]; 

eps = 1e-3; 
maxIter = 2e2; % ����������100�Σ� 
% ����frankworf�㷨 
[sol, x_iters] = frankwolfe(x0,func,grad_func,A,b,Aeq,beq,lb,ub,eps,maxIter,false); 
disp("��������: "+ num2str(size(x_iters,1))); 
disp("�ҵ������Ž�Ϊ:"); 
disp(num2str(sol)); 

% ����ÿһ�ε�����ý�ĺ���ֵ 
y_iters = zeros(1,size(x_iters,1)); 
for i = 1:length(y_iters) 
    y_iters(i) = func(x_iters(i,:)); 
end
figure();
plot(1:length(y_iters), y_iters); 
xlabel('��������');
ylabel('Ŀ�꺯��ֵ'); 

function y = objective(x) 
y = x(1)^2 + x(2)^2 -x(1)*x(2)-2*x(1)+3*x(2); 
end

function gx = grad_objective(x) 
gx = [2*x(1)-x(2)-2; 2*x(2)-x(1)+3; 0; 0]; 
end