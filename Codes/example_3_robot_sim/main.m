clc;
clear;
close all;

tic;
%% Time sequence
T0 = 0;
Ts = 0.001;
Tf = 10;

t = T0:Ts:Tf;

nSteps = numel(t);

%% Controller parameters
eta = 1;
mu = 1;

epsilon = 1e-5;

L1 = 1.5*diag([1,1]);
L2 = 0.2*diag([1,1]);
L_ESO = [L1;L2];

lambda_1 = 10*diag([8,8]);
lambda_2 = 10*diag([9,7]);
lambda_3 = 2*diag([1,5]);

lambda_sw = 0.0001*diag([1,1]);

alpha = 3/5;

h = Ts;
L = 10;
beta = 0.2;

nu = diag([0.003,0.005]);
rho_inf = diag(diag([0.02,0.02]));

m_u_1 = 1*diag([1,1]);

Lu = 2;
Ly = 1;

%% Initialization
nInputs = 2;
nOutputs = 2;

I_n = eye(nOutputs);
I_2n = eye(2*nOutputs);
O_n = zeros(nOutputs);
O_2n = zeros(2*nOutputs);

y_d = zeros(nOutputs,nSteps+1);
y_d(:,1) = [ 0.5*sin((pi/4)*t(1))
             0.2*cos((2*pi/9)*t(1)) ];
y_d(:,2) = [ 0.5*sin((pi/4)*t(2))
             0.2*cos((2*pi/9)*t(2)) ];

q = zeros(nOutputs,nSteps+1);
q(:,1) = [-0.97 1.94]' + 0.01;
q(:,2) = [-0.97 1.94]' + 0.01;

qdot = zeros(nOutputs,nSteps+1);
qdot(:,1) = 0.4;
qdot(:,2) = 0.2;

u = zeros(nInputs,nSteps);
delta_u_eq = zeros(nInputs,nSteps);
delta_u_sw = zeros(nInputs,nSteps);
delta_u = zeros(nInputs,nSteps);
delta_U = zeros(Lu*nInputs,nSteps);

y = zeros(nOutputs,nSteps+1);
y(:,1) = qdot(:,1);
y(:,2) = qdot(:,2);
delta_y = zeros(nOutputs,nSteps+1);
delta_Y = zeros(Ly*nOutputs,nSteps+1);

delta_H = zeros((Ly+Lu)*nOutputs,nSteps);

x_hat = zeros(2*nOutputs,nSteps+1);
x_hat_1 = zeros(nOutputs,nSteps+1);
x_hat_2 = zeros(nOutputs,nSteps+1);

y_hat = zeros(nOutputs,nSteps+1);

e = zeros(nOutputs,nSteps);
e(:,1) = y_d(:,1) - y(:,1);

rho = zeros(nOutputs,nSteps+1);
rho(:,1) = 0.7*[1 1]';
rho(:,2) = 0.7*[1 1]';

tau = zeros(nOutputs,nSteps);

FO_term = zeros(nOutputs,nSteps);
memory_1 = zeros(1,L+1);
memory_1(1) = sig_func(e(1,1),alpha);
memory_2 = zeros(1,L+1);
memory_2(1) = sig_func(e(1,2),alpha);

s = zeros(nOutputs,nSteps);
delta_s = zeros(nOutputs,nSteps);

dist = zeros(nOutputs,nSteps);

PHI_hat = cell(1,nSteps);
PHI_hat{1} = [repmat(diag([10,10]),1,Ly) repmat(diag([1,1]),1,Lu)];
PHI_hat_y = cell(1,nSteps);
PHI_hat_y{1} = repmat(diag([10,10]),1,Ly);
PHI_hat_y_1 = cell(1,nSteps);
PHI_hat_y_1{1} = diag([10,10]);
PHI_hat_u = cell(1,nSteps);
PHI_hat_u{1} = repmat(diag([1,1]),1,Lu);
PHI_hat_u_1 = cell(1,nSteps);
PHI_hat_u_1{1} = diag([1,1]);
PHI_hat_u_others = cell(1,nSteps);
PHI_hat_u_others{1} = repmat(diag([1,1]),1,Lu-1);

fprintf('Simulating PP-MFAFOFTSMC on 2-DOF robotic manipulator...\n');
fprintf('----------------------------------\n');

%% Simulation loop
for k = 2:nSteps
    
    %% PPJM estimation
    PHI_hat{k} = PHI_hat{k-1} + eta*((delta_y(:,k)-(PHI_hat{k-1}+[zeros(nOutputs,Ly*nOutputs) m_u_1 zeros(nOutputs,(Lu-1)*nOutputs)])*delta_H(:,k-1)-x_hat_2(:,k-1))*delta_H(:,k-1)')/(mu+(norm(delta_H(:,k-1),2)^2));
    
    PHI_hat_y{k} = PHI_hat{k}(:,1:Ly*nOutputs);
    PHI_hat_u{k} = PHI_hat{k}(:,Ly*nOutputs+1:end);
    
    PHI_hat_u_1{k} = PHI_hat_u{k}(:,1:nInputs);
    PHI_hat_y_1{k} = PHI_hat_y{k}(:,1:nOutputs);
    
    if(Lu >= 2)
        
        PHI_hat_u_others{k} = PHI_hat_u{k}(:,nInputs+1:end);
        
    end
    
    for ii = 1:nInputs
        for jj = 1:nOutputs
           if(sign(PHI_hat_u_1{k}(ii,jj)) ~= sign(PHI_hat_u_1{1}(ii,jj)))
                PHI_hat_u_1{k}(ii,jj) = PHI_hat_u_1{1}(ii,jj);
           end
        end
    end
    
    if(Lu >= 2)
        
        PHI_hat_u{k} = [PHI_hat_u_1{k} PHI_hat_u_others{k}];
        
    else
        
        PHI_hat_u{k} = PHI_hat_u_1{k};
        
    end
    
    PHI_hat{k} = [PHI_hat_y{k} PHI_hat_u{k}];
    
    %% Reference signals
    y_d(:,k+1) = [ 0.5*sin((pi/4)*t(k))
                   0.2*cos((2*pi/6)*t(k)) ];

    %% Tracking errors
    e(:,k) = y_d(:,k) - y(:,k);
    
    %% Prescribed term
    rho(:,k+1) = (I_n-nu)*rho(:,k) + nu*rho_inf;
    
    %% Transformed tracking errors
    tau(:,k) = (1/2)*log(abs((rho(:,k)+e(:,k))./(rho(:,k)-e(:,k))));
    
    %% Sliding surfaces
    tau_zeta = sig_func(tau(:,k),alpha);
    memory_1 = [tau_zeta(1) memory_1(1:end-1)];
    memory_2 = [tau_zeta(2) memory_2(1:end-1)];
    
    FO_term(1,k) = discrete_FO_coeff_func(beta-1,L,h)*memory_1';
    FO_term(2,k) = discrete_FO_coeff_func(beta-1,L,h)*memory_2';
    
%     s(:,k) = s(:,k-1) + (tau(:,k) - tau(:,k-1)) + lambda_1*tau(:,k) + ...
%              lambda_2*FO_term(:,k-1);
         
    s(:,k) = s(:,k-1) + lambda_3*(tau(:,k) - tau(:,k-1)) + lambda_1*tau(:,k) + ...
             lambda_2*FO_term(:,k-1);
         
    delta_s(:,k) = s(:,k) - s(:,k-1);
    
    %% Temporary variable
%     sigma = exp((((lambda_1+I_n)/2)^-1)*(tau(:,k)-lambda_2*FO_term(:,k)));
    sigma = exp((((lambda_1+lambda_3)/2)^-1)*(lambda_3*tau(:,k)-lambda_2*FO_term(:,k)));
    
    XX{k} = diag(sigma);
    teta_2{k} = (XX{k}^-1)*(XX{k}+eye(2));
    teta_3{k} = XX{k}+eye(2);
    
    %% Controller
    if(Lu >= 2)
        
        delta_u_eq(:,k) = ((PHI_hat_u_1{k}+m_u_1)^-1)*(y_d(:,k+1)-y(:,k)- ...
                           PHI_hat_y{k}*delta_Y(:,k)-PHI_hat_u_others{k}*delta_U(nInputs+1:end,k-1)- ...
                           x_hat_2(:,k)-((sigma-ones(nOutputs,1))./(ones(nOutputs,1)+sigma)).*rho(:,k+1));
    
    else
        
        delta_u_eq(:,k) = ((PHI_hat_u_1{k}+m_u_1)^-1)*(y_d(:,k+1)-y(:,k)- ...
                           PHI_hat_y{k}*delta_Y(:,k)-zeros(nOutputs,1)- ...
                           x_hat_2(:,k)-((sigma-1)./(1+sigma)).*rho(:,k+1));
                   
    end
    
    delta_u_sw(:,k) = ((PHI_hat_u_1{k}+m_u_1)^-1)*(lambda_sw*sign(s(:,k)));
%     delta_u_sw(:,k) = ((PHI_hat_u_1{k}+m_u_1)^-1)*(lambda_sw*sat_func(s(:,k),10));
    
    delta_u(:,k) = delta_u_eq(:,k) + delta_u_sw(:,k);
    
    u(:,k) = u(:,k-1) + delta_u(:,k);
    
    %% Plant simulation
%     dist(:,k) = [ 0.5*sin(qdot(1,k))
%                   0.5*sin(qdot(2,k)) ];
    
    dist(:,k) = [  0.1*sin(t(k)) + 0.003*sin(50*pi*t(k))
                  0.1*cos(2*t(k)) + 0.002*sin(40*pi*t(k)) ];
              
    qddot = twoLink_robot_dynamics(1*t(k),q(:,k),qdot(:,k),u(:,k),dist(:,k));
    qdot(:,k+1) = qdot(:,k) + Ts*qddot;
    q(:,k+1) = q(:,k) + Ts*qdot(:,k+1);
    
    y(:,k+1) = qdot(:,k+1) + dist(:,k);
        
    %% Incremental I/O data
    delta_y(:,k+1) = y(:,k+1) - y(:,k);
    
    delta_U(:,k) = [delta_u(:,k);delta_U(1:(Lu-1)*nInputs,k-1)];
    delta_Y(:,k+1) = [delta_y(:,k+1);delta_Y(1:(Ly-1)*nOutputs,k)];
    
    delta_H(:,k) = [delta_Y(:,k);delta_U(:,k)];
    
    %% ESO
    A = [ I_n  I_n
          O_n  I_n ];
    
    B = [          PHI_hat{k}
          zeros(nOutputs,nOutputs*(Ly+Lu)) ];
    
    C = [ I_n  O_n ];
    
    if(any(abs(eig(A-L_ESO*C)) >= 1))

        disp('Error! ESO matrix is not Hurwitz! Change the parameters of L_ESO.');
        break;
        
    end
    
    x_hat(:,k+1) = A*x_hat(:,k) + B*delta_H(:,k) + L_ESO*(y(:,k)-y_hat(:,k));
    y_hat(:,k+1) = C*x_hat(:,k+1);
    
    x_hat_1(:,k+1) = x_hat(1:nOutputs,k+1);
    x_hat_2(:,k+1) = x_hat(nOutputs+1:end,k+1);
    
end

PHI_hat_u_1_mat = cell2mat(PHI_hat_u_1);
PHI_hat_u_1_mat_11 = PHI_hat_u_1_mat(1,1:nInputs:end);
PHI_hat_u_1_mat_22 = PHI_hat_u_1_mat(2,2:nInputs:end);

PHI_hat_y_1_mat = cell2mat(PHI_hat_y_1);
PHI_hat_y_1_mat_11 = PHI_hat_y_1_mat(1,1:nOutputs:end);
PHI_hat_y_1_mat_22 = PHI_hat_y_1_mat(2,2:nOutputs:end);

y_d = y_d(:,1:end-1);

y = y(:,1:end-1);

rho = rho(:,1:end-1);

fprintf('Successfully done!\n');
fprintf('----------------------------------\n');

%% Analysis
teta2 = cell2mat(teta_2);
teta2_1 = teta2(1,1:2:end);
teta2_2 = teta2(2,2:2:end);

teta3 = cell2mat(teta_3);
teta3_1 = teta3(1,1:2:end);
teta3_2 = teta3(2,2:2:end);

ub_lambda_s_1 = min(min(rho_inf(1,1)*(teta2_1.^-1)),min(rho_inf(1,1)*(teta3_1.^-1)));
ub_lambda_s_2 = min(min(rho_inf(2,1)*(teta2_2.^-1)),min(rho_inf(2,1)*(teta3_2.^-1)));

toc;
%% Plot results
LineWidth_d = 1.5;
LineWidth = 1;

%% PPJM estimation
figure(1);
subplot(1,2,1);
plot(t,PHI_hat_u_1_mat_11,'b','LineWidth',LineWidth);
hold on;
plot(t,PHI_hat_u_1_mat_22,'r--','LineWidth',LineWidth);
grid on;
xlabel('Time (s)','FontName','Times','FontSize',14);
ylabel('$\hat{\Phi}_{u_1}$ (k)','Interpreter','LATEX','FontName','Times','FontSize',14);
% legend('$\hat{\Phi}_{u_{1_{11}}}$','$\hat{\Phi}_{1_{22}}$','Interpreter','LATEX','FontName','Times','FontSize',14);
lg_handle_1_1 = legend('$\hat{\Phi}_{u_{1_{11}}}$','$\hat{\Phi}_{u_{1_{22}}}$');
lg_handle_1_1.Interpreter = 'LATEX';
lg_handle_1_1.FontName = 'Times';
lg_handle_1_1.FontSize = 14;
legend('boxoff');

title('PPJM estimation');

subplot(1,2,2);
plot(t,PHI_hat_y_1_mat_11,'b','LineWidth',LineWidth);
hold on;
plot(t,PHI_hat_y_1_mat_22,'r--','LineWidth',LineWidth);
grid on;
xlabel('Time (s)','FontName','Times','FontSize',14);
ylabel('$\hat{\Phi}_{y_1}$ (k)','Interpreter','LATEX','FontName','Times','FontSize',14);
% legend('$\hat{\Phi}_{y_{1_{11}}}$','$\hat{\Phi}_{y_{1_{22}}}$','Interpreter','LATEX','FontName','Times','FontSize',14);
lg_handle_1_2 = legend('$\hat{\Phi}_{y_{1_{11}}}$','$\hat{\Phi}_{y_{1_{22}}}$');
lg_handle_1_2.Interpreter = 'LATEX';
lg_handle_1_2.FontName = 'Times';
lg_handle_1_2.FontSize = 14;
legend('boxoff');

%% Difference of sliding surfaces
figure(2);
subplot(2,1,1);
plot(t,delta_s(1,:),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (s)','FontName','Times','FontSize',14);
ylabel('$\Delta{s}_1$ (rad/s)','Interpreter','LATEX','FontName','Times','FontSize',14);

title('Difference of sliding surfaces');

subplot(2,1,2);
plot(t,delta_s(2,:),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (s)','FontName','Times','FontSize',14);
ylabel('$\Delta{s}_2$ (rad/s)','Interpreter','LATEX','FontName','Times','FontSize',14);

%% Tracking errors
figure(3);
subplot(2,1,1);
plot(t,rho(1,:),'m--','LineWidth',LineWidth_d);
hold on;
plot(t,e(1,:),'b','LineWidth',LineWidth);
plot(t,-rho(1,:),'m--','LineWidth',LineWidth_d);
grid on;
xlabel('Time (s)','FontName','Times','FontSize',14);
ylabel('$e_1$ (rad/s)','Interpreter','LATEX','FontName','Times','FontSize',14);
% legend('Upper prescribed boundary','Tracking error','Lower prescribed boundary','FontName','Times','FontSize',14);
lg_handle_2_1 = legend('Upper prescribed boundary','Tracking error','Lower prescribed boundary');
lg_handle_2_1.FontName = 'Times';
lg_handle_2_1.FontSize = 14;
axes('Position',[0.271595900439239 0.60989247311828 0.163250366032211 0.098310291858679]);
interval_1 = 1000:3500;
plot(t(interval_1),rho(1,interval_1),'m--','LineWidth',LineWidth_d);
hold on;
plot(t(interval_1),e(1,interval_1),'b','LineWidth',LineWidth);
plot(t(interval_1),-rho(1,interval_1),'m--','LineWidth',LineWidth_d);
xlim([t(interval_1(1)) t(interval_1(end))]);
ylim([-rho_inf(1)-rho_inf(1)/1 rho_inf(1)+rho_inf(1)/1]);

title('Tracking errors');

subplot(2,1,2);
plot(t,rho(2,:),'m--','LineWidth',LineWidth_d);
hold on;
plot(t,e(2,:),'b','LineWidth',LineWidth);
plot(t,-rho(2,:),'m--','LineWidth',LineWidth_d);
grid on;
xlabel('Time (s)','FontName','Times','FontSize',14);
ylabel('$e_2$ (rad/s)','Interpreter','LATEX','FontName','Times','FontSize',14);
% legend('Upper prescribed boundary','Tracking error','Lower prescribed boundary','FontName','Times','FontSize',14);
lg_handle_2_2 = legend('Upper prescribed boundary','Tracking error','Lower prescribed boundary');
lg_handle_2_2.FontName = 'Times';
lg_handle_2_2.FontSize = 14;
axes('Position',[0.27086383601757 0.135544487799413 0.163250366032211 0.085500592787076]);
interval_2 = 1000:3500;
plot(t(interval_2),rho(2,interval_1),'m--','LineWidth',LineWidth_d);
hold on;
plot(t(interval_2),e(2,interval_2),'b','LineWidth',LineWidth);
plot(t(interval_2),-rho(2,interval_2),'m--','LineWidth',LineWidth_d);
xlim([t(interval_2(1)) t(interval_2(end))]);
ylim([-rho_inf(2)-rho_inf(2)/1 rho_inf(2)+rho_inf(2)/1]);

%% Control inputs
figure(4);
subplot(2,1,1);
plot(t,u(1,:),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (s)','FontName','Times','FontSize',14);
ylabel('$u_1$ (N.m)','Interpreter','LATEX','FontName','Times','FontSize',14);

title('Control inputs');

subplot(2,1,2);
plot(t,u(2,:),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (s)','FontName','Times','FontSize',14);
ylabel('$u_2$ (N.m)','Interpreter','LATEX','FontName','Times','FontSize',14);

%% Tracking performance
figure(5);
subplot(2,1,1);
plot(t,y_d(1,:),'r--','LineWidth',LineWidth_d);
hold on;
plot(t,y(1,:),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (s)','FontName','Times','FontSize',14);
ylabel('$\dot{q}_1$ (rad/s)','Interpreter','LATEX','FontName','Times','FontSize',14);
% legend('Desired','Actual','Location','SouthEast','FontName','Times','FontSize',14);
lg_handle_3_1 = legend('Desired','Actual');
lg_handle_3_1.Location = 'SouthEast';
lg_handle_3_1.FontName = 'Times';
lg_handle_3_1.FontSize = 14;
legend('boxoff');
title('Angular velocities','FontName','Times','FontSize',14);

title('Tracking performance');

subplot(2,1,2);
plot(t,y_d(2,:),'r--','LineWidth',LineWidth_d);
hold on;
plot(t,y(2,:),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (s)','FontName','Times','FontSize',14);
ylabel('$\dot{q}_2$ (rad/s)','Interpreter','LATEX','FontName','Times','FontSize',14);
% legend('Desired','Actual','Location','SouthEast','FontName','Times','FontSize',14);
lg_handle_3_2 = legend('Desired','Actual');
lg_handle_3_2.Location = 'NorthEast';
lg_handle_3_2.FontName = 'Times';
lg_handle_3_2.FontSize = 14;
legend('boxoff');
