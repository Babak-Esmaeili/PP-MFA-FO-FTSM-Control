clc;
clear;
close all;

tic;
%% Time sequence
T0 = 0;
Ts = 0.001;
Tf = 10;

% t = T0:Ts:Tf;
t = 1:500;

% nSteps = numel(t);
nSteps = 500;

%% Controller parameters
eta = 0.8;
mu = 0.5;

epsilon = 1e-5;

L1 = 0.5*1;
L2 = 0.2*1;
L_ESO = [L1;L2];

lambda_1 = 20*1;
lambda_2 = 0.1*1;
lambda_3 = 04*1;

lambda_sw = 0.0015*1;

alpha = 3/5;

h = Ts;
L = 8;
beta = 0.4;

nu = 0.1;
rho_inf = 0.03;

m_u_1 = 1*1;

Lu = 1;
Ly = 1;

%% Initialization
nInputs = 1;
nOutputs = 1;

I_n = eye(nOutputs);
I_2n = eye(2*nOutputs);
O_n = zeros(nOutputs);
O_2n = zeros(2*nOutputs);

u = zeros(nInputs,nSteps);
delta_u_eq = zeros(nInputs,nSteps);
delta_u_sw = zeros(nInputs,nSteps);
delta_u = zeros(nInputs,nSteps);
delta_U = zeros(Lu*nInputs,nSteps);

y = zeros(nOutputs,nSteps+1);
y(:,1) = 0;
y(:,2) = 0;
delta_y = zeros(nOutputs,nSteps+1);
delta_Y = zeros(Ly*nOutputs,nSteps+1);

delta_H = zeros((Ly+Lu)*nOutputs,nSteps);

x_hat = zeros(2*nOutputs,nSteps+1);
x_hat_1 = zeros(nOutputs,nSteps+1);
x_hat_2 = zeros(nOutputs,nSteps+1);

y_hat = zeros(nOutputs,nSteps+1);

y_d = zeros(nOutputs,nSteps+1);
y_d(1) = 0.2*(sin((2*pi/50)*(0+1))+sin((2*pi/100)*(0+1))+sin((2*pi/150)*(0+1)));
y_d(2) = 0.2*(sin((2*pi/50)*(1+1))+sin((2*pi/100)*(1+1))+sin((2*pi/150)*(1+1)));
e = zeros(nOutputs,nSteps);
e(:,1) = y_d(:,1) - y(:,1);

rho = zeros(nOutputs,nSteps+1);
rho(:,1) = 1.2;
rho(:,2) = 1.2;

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
PHI_hat{1} = [repmat(0.5*1,1,Ly) repmat(0.5*1,1,Lu)];
PHI_hat_y = cell(1,nSteps);
PHI_hat_y{1} = repmat(0.5*1,1,Ly);
PHI_hat_y_1 = cell(1,nSteps);
PHI_hat_y_1{1} = 0.5*1;
PHI_hat_u = cell(1,nSteps);
PHI_hat_u{1} = repmat(0.5*1,1,Lu);
PHI_hat_u_1 = cell(1,nSteps);
PHI_hat_u_1{1} = 0.5*1;
PHI_hat_u_others = cell(1,nSteps);
PHI_hat_u_others{1} = repmat(0.5*1,1,Lu-1);

fprintf('Simulating PP-MFAFOFTSMC on a numerical example...\n');
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
    y_d(:,k+1) = 0.2*(sin((2*pi/50)*(k+1))+sin((2*pi/100)*(k+1))+sin((2*pi/150)*(k+1)));
    
    %% Tracking errors
    e(:,k) = y_d(:,k) - y(:,k);
    
    %% Prescribed term
    rho(:,k+1) = (I_n-nu)*rho(:,k) + nu*rho_inf;
    
    %% Transformed tracking errors
    tau(:,k) = (1/2)*log(abs((rho(:,k)+e(:,k))./(rho(:,k)-e(:,k))));
    
    %% Sliding surfaces
    tau_zeta = sig_func(tau(:,k),alpha);
    memory_1 = [tau_zeta(1) memory_1(1:end-1)];
    
    FO_term(1,k) = discrete_FO_coeff_func(beta-1,L,h)*memory_1';
    
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
    dist(:,k) = 0.5 + 0.15*sin(k/30);
              
    y(:,k+1) = (y(k)/(1+y(k)^2)) + (u(k)^3) + dist(:,k);
    
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

PHI_hat_y_1_mat = cell2mat(PHI_hat_y_1);
PHI_hat_y_1_mat_11 = PHI_hat_y_1_mat(1,1:nInputs:end);

y_d = y_d(:,1:end-1);

y = y(:,1:end-1);

rho = rho(:,1:end-1);

fprintf('Successfully done!\n');
fprintf('----------------------------------\n');

%% Analysis
teta2 = cell2mat(teta_2);
teta2_1 = teta2(1,1:2:end);

teta3 = cell2mat(teta_3);
teta3_1 = teta3(1,1:2:end);

ub_lambda_s_1 = min(min(rho_inf(1,1)*(teta2_1.^-1)),min(rho_inf(1,1)*(teta3_1.^-1)));

%% Indices
IAE = sum(abs(e'))'*1
ISE = sum(abs(e').^2)'*1
ITAE = sum(repmat(t',1,nOutputs).*abs(e'))'
TV = sum(abs(delta_u'))'*1
L2 = sqrt(sum(u'.^2))*1

toc;
%% Plot results
LineWidth_d = 1.5;
LineWidth = 1;

%% PPJM estimation
% figure(1);
% subplot(1,2,1);
% plot(t,PHI_hat_u_1_mat_11,'b','LineWidth',LineWidth);
% hold on;
% plot(t,PHI_hat_u_1_mat_22,'r--','LineWidth',LineWidth);
% grid on;
% xlabel('Time (s)','FontName','Times','FontSize',14);
% ylabel('$\hat{\Phi}_{u_1}$ (k)','Interpreter','LATEX','FontName','Times','FontSize',14);
% % legend('$\hat{\Phi}_{u_{1_{11}}}$','$\hat{\Phi}_{1_{22}}$','Interpreter','LATEX','FontName','Times','FontSize',14);
% lg_handle_1_1 = legend('$\hat{\Phi}_{u_{1_{11}}}$','$\hat{\Phi}_{u_{1_{22}}}$');
% lg_handle_1_1.Interpreter = 'LATEX';
% lg_handle_1_1.FontName = 'Times';
% lg_handle_1_1.FontSize = 14;
% legend('boxoff');
% 
% subplot(1,2,2);
% plot(t,PHI_hat_y_1_mat_11,'b','LineWidth',LineWidth);
% hold on;
% plot(t,PHI_hat_y_1_mat_22,'r--','LineWidth',LineWidth);
% grid on;
% xlabel('Time (s)','FontName','Times','FontSize',14);
% ylabel('$\hat{\Phi}_{y_1}$ (k)','Interpreter','LATEX','FontName','Times','FontSize',14);
% % legend('$\hat{\Phi}_{y_{1_{11}}}$','$\hat{\Phi}_{y_{1_{22}}}$','Interpreter','LATEX','FontName','Times','FontSize',14);
% lg_handle_1_2 = legend('$\hat{\Phi}_{y_{1_{11}}}$','$\hat{\Phi}_{y_{1_{22}}}$');
% lg_handle_1_2.Interpreter = 'LATEX';
% lg_handle_1_2.FontName = 'Times';
% lg_handle_1_2.FontSize = 14;
% legend('boxoff');
% 
% suptitle('PPJM estimation');

%% Difference of sliding surface
figure(2);
plot(t,delta_s(1,:),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (s)','FontName','Times','FontSize',14);
ylabel('$\Delta{s}_1$ (k)','Interpreter','LATEX','FontName','Times','FontSize',14);

title('Difference of sliding surface');

%% Tracking error
figure(3);
plot(t,rho(1,:),'m--','LineWidth',LineWidth_d);
hold on;
plot(t,-e(1,:),'b','LineWidth',LineWidth);
plot(t,-rho(1,:),'m--','LineWidth',LineWidth_d);
grid on;
xlabel('Time (s)','FontName','Times','FontSize',14);
ylabel('$e_1$ (k)','Interpreter','LATEX','FontName','Times','FontSize',14);
% legend('Upper prescribed boundary','Tracking error','Lower prescribed boundary','FontName','Times','FontSize',14);
lg_handle_2_1 = legend('Upper prescribed boundary','Tracking error','Lower prescribed boundary');
lg_handle_2_1.FontName = 'Times';
lg_handle_2_1.FontSize = 14;
axes('Position',[0.271595900439239 0.60989247311828 0.163250366032211 0.098310291858679]);
interval_1 = 40:60;
plot(t(interval_1),rho(1,interval_1),'m--','LineWidth',LineWidth_d);
hold on;
plot(t(interval_1),e(1,interval_1),'b','LineWidth',LineWidth);
plot(t(interval_1),-rho(1,interval_1),'m--','LineWidth',LineWidth_d);
xlim([t(interval_1(1)) t(interval_1(end))]);
ylim([-rho_inf(1)-rho_inf(1)/1 rho_inf(1)+rho_inf(1)/1]);

title('Tracking error');

%% Control input
figure(4);
plot(t,u(1,:),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (s)','FontName','Times','FontSize',14);
ylabel('$u_1$ (k)','Interpreter','LATEX','FontName','Times','FontSize',14);

title('Control input');

%% Tracking performance
figure(5);
plot(t,y_d(1,:),'r--','LineWidth',LineWidth_d);
hold on;
plot(t,y(1,:),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (s)','FontName','Times','FontSize',14);
ylabel('$y_1$ (k)','Interpreter','LATEX','FontName','Times','FontSize',14);
% legend('Desired','Actual','Location','SouthEast','FontName','Times','FontSize',14);
lg_handle_3_1 = legend('Desired','Actual');
lg_handle_3_1.Location = 'SouthEast';
lg_handle_3_1.FontName = 'Times';
lg_handle_3_1.FontSize = 14;
legend('boxoff');

title('Tracking performance');
