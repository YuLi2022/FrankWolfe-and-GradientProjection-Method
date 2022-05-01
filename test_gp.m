% % 梯度投影法：计算示例 
% func = @objective; 
% grad_func = @grad_objective; 
% A = [-1,0; 0,-1; 1,0; 0,1]; 
% b = [-5; -8; 0; 0]; 
% Aeq = []; 
% beq = [];  
% x0=[0;0]; 
% eps = 1e-6; 
% maxIter = 1e3; 
% 
% % 调用梯度投影法 
% [sol, x_iters] = GradientProjection(x0, func, grad_func, A, b, Aeq, beq, eps, maxIter, true); 
% disp("==================="); 
% disp("迭代次数: "+ num2str(size(x_iters,1))); 
% disp("找到的最优解为:"); 
% disp(num2str(sol)); 
% 
% % 绘图 
% y_iters = zeros(1,size(x_iters,1)); 
% for i = 1:length(y_iters) 
%     y_iters(i) = func(x_iters(i,:)); 
% end
% figure();
% plot(y_iters); 
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


%测试函数1： 
func = @objective; 
grad_func = @grad_objective; 
A = [-1,-1;-1,-5;1,0;0,1]; 
b = [-2;-5;0;0]; 
Aeq = []; 
beq = []; 
x0=[0;0];
eps = 1e-6; 
maxIter = 1e3; 

% 调用梯度投影法 
[sol, x_iters] = GradientProjection(x0, func, grad_func, A, b, Aeq, beq, eps, maxIter, true); 
disp("==================="); 
disp("迭代次数: "+ num2str(size(x_iters,1))); 
disp("找到的最优解为:"); 
disp(num2str(sol)); 

% 绘图 
y_iters = zeros(1,size(x_iters,1)); 
for i = 1:length(y_iters) 
    y_iters(i) = func(x_iters(i,:)); 
end
figure();
plot(y_iters); 
xlabel('迭代次数'); 
ylabel('目标函数值'); 

function y = objective(x) 
y = 2*x(1)^2 + 2*x(2)^2 - 2*x(1)*x(2) - 4*x(1) - 6*x(2); 
end

function gx = grad_objective(x) 
gx = [4*x(1)-2*x(2)-4; 4*x(2)-2*x(1)-6]; 
end

% 测试函数2 
% func = @objective; 
% grad_func = @grad_objective; 
% A = [1,-2,0; 1,0,0; 0,1,0; 0,0,1]; 
% b = [-3;0;0;0]; 
% Aeq = [1,1,1]; 
% beq = [2]; 
% x0=[1;0;1]; 
% eps = 1e-5; 
%  
% sol = GradientProjection(x0, func, grad_func, A, b, Aeq, beq, eps) ; 
% 
% function y = objective(x) 
% y = x(1)^2 + x(1)*x(2) + 2*x(2)^2 - 6*x(1) - 2*x(2) - 12*x(3); 
% end
% 
% function gx = grad_objective(x) 
% gx = [2*x(1)+x(2)-6; x(1)+4*x(2)-2; -12]; 
% end