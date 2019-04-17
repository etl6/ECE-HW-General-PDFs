%% Emma LaPorte (etl6)
% ECE 681 
% HW1 - ROC Generating Problems

%% for a new dataset need to change the filename that goes into csvread() 
 filename = 'logNormalData.csv';
 
 
 table = csvread(filename);
 tableSorted = sortrows(table,2); % sorts ascending according to lambda
 truthSorted = tableSorted(:,1); % col vector
 lambdaSorted = tableSorted(:,2); % col vector
 lambda_H0 = lambdaSorted(truthSorted==0); % H0 decision statistics
 lambda_H1 = lambdaSorted(truthSorted==1); % H1 decision statistics

% case 1
 Beta1=[-inf; lambdaSorted; inf];
 [Pf_1,Pd_1] = calcPfPd(Beta1, lambda_H0, lambda_H1);
 figure()
 plot(Pf_1,Pd_1,'linewidth',5)
 hold on
% case 2
 Beta2=linspace(min(lambdaSorted),max(lambdaSorted),99)';
 Beta2=[-Inf; Beta2; Inf];
 [Pf_2,Pd_2] = calcPfPd(Beta2, lambda_H0, lambda_H1); 
 plot(Pf_2,Pd_2)
% case 3 
 if length(lambdaSorted) < 99
     Beta3 = [-Inf; lambdaSorted; Inf];
 else 
    index3 = (round(length(lambdaSorted)*(1:99)/99))';
    Beta3 = [-Inf; lambdaSorted(index3); Inf];
 end
 [Pf_3,Pd_3] = calcPfPd(Beta3, lambda_H0, lambda_H1); 
 plot(Pf_3,Pd_3,'--')
% case 4
 Beta4 = [-Inf; lambda_H0; Inf];
 [Pf_4,Pd_4] = calcPfPd(Beta4, lambda_H0, lambda_H1); 
 plot(Pf_4,Pd_4,':')
% case 5
 if length(lambda_H0)>101
    increment = round(length(lambda_H0)/101);
    Beta5 = [-Inf; lambda_H0(1:increment:end); Inf];
 else
     Beta5 = [-Inf; lambda_H0; Inf]; % for smaller datasets
 end
 [Pf_5,Pd_5] = calcPfPd(Beta5, lambda_H0, lambda_H1); 
 plot(Pf_5,Pd_5)

 % testing perfcurve MATLAB function
 [X,Y,BetaPerfCurve] = perfcurve(truthSorted,lambdaSorted,1);
 plot(X,Y)
 
 hold off
 legend('Case 1','Case 2','Case 3','Case 4','Case 5','perfcurve MATLAB')
 title('ROCs (etl6)')

 %% Emma LaPorte
% ECE 681
% Problem 12

table = csvread('rocData.csv');
Pf = table(:,1); % col vector
Pd = table(:,2); % col vector
figure()
plot(Pf,Pd)
hold on
title('ROC and Operating Points for rocData.csv (etl6)')
ylabel('P_{D}');xlabel('P_{FA}');

Pcd = 0.5*Pd + 0.5*(1-Pf);
[Pcd1,index1] = max(Pcd);
plot(Pf(index1),Pd(index1),'r*')
Pcd_str = num2str(Pcd1);

Pcd = (2/3)*Pd + (1/3)*(1-Pf);
[Pcd2,index2] = max(Pcd);
plot(Pf(index2),Pd(index2),'k*')
Pcd_str2 = num2str(Pcd2);

Pcd = (1/3)*Pd + (2/3)*(1-Pf);
[Pcd3,index3] = max(Pcd);
plot(Pf(index3),Pd(index3),'g*')
Pcd_str3 = num2str(Pcd3);

AUC = trapz(Pf,Pd);
legend('ROC curve','max P_{cd} when p(H_{0})=p(H_{1})','max P_{cd} when 2p(H_{0})=p(H_{1})','max P_{cd} when p(H_{0})=2p(H_{1})')

hold off
 

%% function for calculating Pf, Pd pairs
function [Pf,Pd] = calcPfPd(Beta1, lambda_H0, lambda_H1)
     Pd = zeros(length(Beta1),1);
     Pf = zeros(length(Beta1),1);
     for i = 1:length(Beta1) % iterating thru Beta
         Pf_num=0;
         Pd_num=0;
         for j=1:length(lambda_H0)
             if lambda_H0(j)>=Beta1(i)
                 Pf_num=Pf_num+1;
             end
         end
         for k=1:length(lambda_H1)
             if lambda_H1(k)>=Beta1(i)
                 Pd_num=Pd_num+1;
             end
         end
         Pd(i) = Pd_num/length(lambda_H1);
         Pf(i) = Pf_num/length(lambda_H0);
     end
end

