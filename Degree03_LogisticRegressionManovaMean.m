clear all
close all
clc

addpath('Resources') 
%% Generate folder for results
folder_name='Results/Degree03_LogisticRegressionManovaMean'; %Name for the results folder: it should be named after the kind of test performed

currDate = datestr(datetime,30);
mkdir(folder_name,currDate);
results_folder=strcat(folder_name,'/',currDate);
%All the figures we want to save should be added to the array of figures
%through figs(1)=figure(1); figs(2)=figure(2);

%% Load the dataset
% Load the dataset from ProjectionAndOthers-MULTIPLE
% [counter r1 r2 r3 gamma Gauss_test_HZ MBox_p Manova_d Manova_p]
load('Results/Degree03_NoiseCircular_SimulManova/20210809T235258/dataset.mat');
gamma = dataset(:,5);
Gauss_test_HZ = dataset(:,6);
MBox_p = dataset(:,7);
Manova_d = dataset(:,8);
Manova_p = dataset(:,9);


%% Setting the model variables
x = log(abs(gamma));
t = (Manova_p >= 0.05);

%% Logistic regression
t = t+1;
[B, ~, stats] = mnrfit(x,t); % computes the weight matrix

lev=0.46; %Model threshold
sep_line=(log((1-lev)./lev)-B(1))/B(2);
exp(sep_line)

pihat = mnrval(B,x);
t_pred = pihat(:,2)>=lev; %Predicted value
t_pred=t_pred+1;

cm=confusionmat(t,t_pred)
fpr=cm(1,2)/(cm(1,2)+cm(1,1)); %1-specificity
tpr=cm(2,2)/(cm(2,1)+cm(2,2)); %sensitivuty


%% Plots
figs(1)=figure(1);
plot(x,Manova_p,'x'); hold on; grid on;
% plot(x,t-1,'rx');
yline(0.05,'r');
plot(x,pihat(:,2),'r.');
xline(sep_line,'b--');
legend("P-value of MANOVA test","Level $\alpha=0.05$ for MANOVA test","Fitted model","Separating value $\overline{\gamma}_{MANOVA}$","Location","Northwest","interpreter","latex");
title("MANOVA test");
xlabel("log(\gamma(z_0))");
ylabel("Probability");
hold off

figs(2)=figure(2);
plot(Manova_p,pihat(:,2),'x'); hold on; grid on;
yline(lev,'b--');
xline(0.05,'r');
xlabel("P-value of MANOVA test");
ylabel("Probability from fitted model");
legend("Fitted vs actual","Chosen model threshold","Level $\alpha=0.05$ for MANOVA test","Location","Southeast","interpreter","latex")
title("Logistic regression fitted probabilities vs MANOVA test P-values");
hold off

figs(3)=figure(3); % ROC curve
[X,Y]=perfcurve(t,x,2);
plot(X,Y);
xlabel('False positive rate') 
ylabel('True positive rate')
title('ROC for Classification by Logistic Regression')
hold on; grid on
plot(fpr,tpr,'rx');


%% Save workspace and figures to the folder
savefig(figs,strcat(results_folder,'/figures.fig'),'compact');
clear figs
save(strcat(results_folder,'/workspace'));