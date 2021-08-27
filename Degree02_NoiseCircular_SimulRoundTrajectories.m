clear all
close all
clc

addpath('Resources') 
%% Generate folder for results
folder_name='Results/Degree02_NoiseCircular_SimulRoundTrajectories'; %Name for the results folder: it should be named after the kind of test performed

currDate = datestr(datetime,30);
mkdir(folder_name,currDate);
results_folder=strcat(folder_name,'/',currDate);
%All the figures we want to save should be added to the array of figures
%through figs(1)=figure(1); figs(2)=figure(2);

%% Parameters
N=2; % order of the polynomial
sigma_a=.01; % variance of the noise on coefficients
K=2; % Number of iterations per simulation (n of noisy measurements per polynomial)
scale=2; % the roots of the polynomial will be generated with Re(r),Im(r) in [-scale +scale], a square
T=50; % Number of poits in each trajectory
lambdas=linspace(0,1,T);

%% Generate Polynomial and Covariance matrix
% Generate a random circular covariance
[Sigma,C_atilda,A] = generate_covariance(N,sigma_a,'circular');
% Generate a polynomial
a=[1; 0; -1];
% Compute corresponding roots
r=roots(a);

%% Simulation
h = waitbar(0,'Simulations in progress ... Please wait...');
a_n=zeros(N+1,K,T); %Matrix to contain coefficients at every iteration
r_n=zeros(N,K,T); %Matrix to collect roots computed at every iteration
a_mod=[1; 0; 0.5+0.5*1i];%[1; 0; 1.5+0.5*1i]; %Add noise to coefficients
r_mod=roots(a_mod);
r_prec=r;
a_n(:,1,:)=curved_trajectory1_vector(a,a_mod,lambdas);
a_n(:,2,:)=curved_trajectory2_vector(a,a_mod,lambdas);
for k=1:K
%     noise_tilda=A*randn(2*N,1); %Generate colored noise
    for t=1:T
        r_curr=roots(a_n(:,k,t));
        r_n(:,k,t)=r_curr(order_roots_permutations(r_curr,r_prec)); % Order the estimated roots
        r_prec=r_n(:,k,t);
        
        waitbar((k-1)*T+t/(K*T)) %Update waitbar
    end
end

close(h); %Close waitbar

%% Plots
figs(1)=figure(1);

subplot(1,3,1);
%viscircles([0 0],1,'color','b','linestyle','--','LineWidth',0.1);
hold on;
plot(zeros(2,1),5*[-1,1],'b--','LineWidth',0.1);plot(5*[-1,1],zeros(2,1),'b--','LineWidth',0.1);
for k=1:K
    for ii=2:N+1
        plot(real(reshape(a_n(ii,k,:),[1,T])),imag(reshape(a_n(ii,k,:),[1,T])),'.','MarkerSize',5); hold on; % Simulated coeff
    end
end
% for ii=1:N
%     ellipse_plot(0.1*inv(C_atilda([ii N+ii],[ii N+ii])),[real(a(1+ii)),imag(a(1+ii))])
% end
plot(real(a),imag(a),'*k','MarkerSize',20);
plot(real(a_mod(3)),imag(a_mod(3)),'*g','MarkerSize',20);
axis equal;axis(1.5*[-1,1,-1,1]);
title("Coefficients");grid on;hold off

subplot(1,3,2);
%viscircles([0 0],1,'color','b','linestyle','--','LineWidth',0.1);
hold on;
plot(zeros(2,1),5*[-1,1],'b--','LineWidth',0.1);plot(5*[-1,1],zeros(2,1),'b--','LineWidth',0.1);
for k=1
    for t=1:T
        plot(real(r_n(1,k,t)),imag(r_n(1,k,t)),'r.','MarkerSize',5); hold on; % Simulated roots
        plot(real(r_n(2,k,t)),imag(r_n(2,k,t)),'b.','MarkerSize',5); hold on;
    end
end
plot(real(r),imag(r),'*k','MarkerSize',20); % True roots
plot(real(r_mod),imag(r_mod),'*g','MarkerSize',20); % Modified roots
axis equal;axis(1.5*[-1,1,-1,1]);
title("Roots");grid on;hold off

subplot(1,3,3);
%viscircles([0 0],1,'color','b','linestyle','--','LineWidth',0.1);
hold on;
plot(zeros(2,1),5*[-1,1],'b--','LineWidth',0.1);plot(5*[-1,1],zeros(2,1),'b--','LineWidth',0.1);
for k=2
    for t=1:T
        plot(real(r_n(1,k,t)),imag(r_n(1,k,t)),'r.','MarkerSize',5); hold on; % Simulated roots
        plot(real(r_n(2,k,t)),imag(r_n(2,k,t)),'b.','MarkerSize',5); hold on;
    end
end
plot(real(r),imag(r),'*k','MarkerSize',20); % True roots
plot(real(r_mod),imag(r_mod),'*g','MarkerSize',20); % Modified roots
axis equal;axis(1.5*[-1,1,-1,1]);
title("Roots");grid on;hold off

%% Save workspace and figures to the folder
savefig(figs,strcat(results_folder,'/figures.fig'),'compact');
clear figs
save(strcat(results_folder,'/workspace'));