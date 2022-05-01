% % frank worf算法测试算例 
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
% maxIter = 2e3; % 最大允许迭代100次； 
% % 调用frankworf算法 
% [sol, x_iters] = frankwolfe(x0,func,grad_func,A,b,Aeq,beq,lb,ub,eps,maxIter,false); 
% disp("迭代次数: "+ num2str(size(x_iters,1))); 
% disp("找到的最优解为:"); 
% disp(num2str(sol)); 
% 
% % 计算每一次迭代获得解的函数值 
% y_iters = zeros(1,size(x_iters,1)); 
% for i = 1:length(y_iters) 
%     y_iters(i) = func(x_iters(i,:)); 
% end
% figure();
% plot(1:length(y_iters), y_iters); 
% xlabel('迭代次数');
% ylabel('目标函数值'); 
% 
% % 定义目标函数和梯度函数 
% function y = objective(x) 
% y = 2*x(1)^2 + 3*x(2)^2 - 4*x(1)*x(2)-10*x(1) + 7; 
% end
% 
% function gx = grad_objective(x) 
% gx = [4*x(1)-4*x(2)-10; 6*x(2)-4*x(1)]; 
% end



func = @objective; 
grad_func = @grad_objective; 

% 测试函数  
A = []; 
b = []; 
Aeq = [1,1,1,0;1,5,0,1]; 
beq = [3;6]; 
lb = [0,0,0,0];  
ub = []; 
%x0 = [2;0;1;4];
x0 = [0.5;1;0.5;0.5]; 

eps = 1e-3; 
maxIter = 2e2; % 最大允许迭代100次； 
% 调用frankworf算法 
[sol, x_iters] = frankwolfe(x0,func,grad_func,A,b,Aeq,beq,lb,ub,eps,maxIter,false); 
disp("迭代次数: "+ num2str(size(x_iters,1))); 
disp("找到的最优解为:"); 
disp(num2str(sol)); 

% 计算每一次迭代获得解的函数值 
y_iters = zeros(1,size(x_iters,1)); 
for i = 1:length(y_iters) 
    y_iters(i) = func(x_iters(i,:)); 
end
figure();
plot(1:length(y_iters), y_iters); 
xlabel('迭代次数');
ylabel('目标函数值'); 

function y = objective(x) 
y = x(1)^2 + x(2)^2 -x(1)*x(2)-2*x(1)+3*x(2); 
end

function gx = grad_objective(x) 
gx = [2*x(1)-x(2)-2; 2*x(2)-x(1)+3; 0; 0]; 
end