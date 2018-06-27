% Script to Reify 
close all
clear all
%% Fit Piecwise Polynomial to Experimental Tests to Select Ground Truth
%Fit and Plot Stress Strain Data Using Cubic Spline

filename = 'StressStrainImportTest';
sheet = ["17-080","17-167", "17-162", "17-166", "17-172"];

for i = 1:length(sheet)
    A = xlsread(filename, sheet(i));
    A = [A; A(end,1)+0.001,0];
    %A = [A; linspace(A(end,1)+0.002,0.14,10)', linspace(A(end,2)-40,0,10)']; %exp(-linspace(A(end,1)+0.002,0.14,10)')]  %zeros(10,1)];    %concatenate zeros to the end of this variable by making
    xq = linspace(0,A(end,1),length(A(:,1)));
    y = A(:,2);
    s(i) = spline(A(:,1), y);
    p = plot(A(:,1),A(:,2),'o',xq,ppval(s(i),xq),'-');
%     lgd(i) = legend(strcat('y = 17-', sheet(1,i)));
    hold on
end

l = linspace(0,0.1226); %0.083); %

for i = 1:length(l)
    for j = 1:length(s)
        py(j,i) = ppval(s(j),l(1,i));
    end
    avg_y(1,i) = mean(py(:,i));
end

% Plot Average Values and Peicewise Polynomials
    figure
    plot(l,avg_y, 'o k')
    hold on
    
    for i = 1:length(s)
     plot(l, ppval(s(i),l))
     hold on
    end 

str = 'Stress \sigma vs. Strain \epsilon';
text(0.1,1800,str,'Interpreter','tex')
title(['Stress \sigma vs. Strain \epsilon'],'Interpreter','tex')

xlabel(['Strain \epsilon'], 'Interpreter','tex')
ylabel(['Stress \sigma'], 'Interpreter','tex')

j=1;
B= cell(1,5);
for i = 1:length(sheet)
    B{i} = xlsread(filename, sheet(i));
end
 
%% Apply a gp to the Data Set to Find Mean and Variance of Each
meanfunc = {'meanZero'};
covfunc = {'covSEiso'}; ell = 0.1; sf=0.1;
likfunc = {'likGauss'}; sn = 1.742;      

hyp = struct('mean', [], 'cov', log([ell; sf]), 'lik', log(sn));

hyp2 = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, B{1,i}(:,1), B{1,i}(:,2));

for i = 1:length(sheet)
    xs = linspace(0, B{1,i}(end,1), 100)';
    [mu{i} s2{i}] = gp(hyp2, @infGaussLik, meanfunc, covfunc, likfunc, B{1,i}(:,1), B{1,i}(:,2), xs);

    f = [mu{1,i}(:,1)+2*sqrt(s2{1,i}(:,1)); flipdim(mu{1,i}(:,1)-2*sqrt(s2{1,i}(:,1)),1)];
    fill([xs; flipdim(xs,1)], f, [7 7 7]/8)
    hold on; plot(xs, mu{1,i}(:,1)); plot(B{1,i}(:,1), B{1,i}(:,2), '+')
end

%% Add Observational Error to Experimental Sources

%Plot of Average of Spline Means
l_mod=l';
xs_test = linspace(0, 0.1226)'; %0.083)';%%max(cellfun(@(c) c(:,1), B(1,:),1:5,'un',0)))' %B{1,1}(end,1), 100)'

r = ppval(s(1),xs_test);

for i = 1:length(l) %length(mu{1,1}(:,1))      %It goes the length of the shortest one 
    med_of_spl(i,1) = median(py(:,i)); % median(cellfun(@(c) c(i,1), mu(1,:))); %  

    mmax(i,1) = max(py(:,i));       %cellfun(@(c) c(i,1), mu(1,:)));
    mmin(i,1) = min(py(:,i));        %cellfun(@(c) c(i,1), mu(1,:)));
    
    mdiff(i,1) = mmax(i,1)- mmin(i,1);
    
    
end

    f2 = [med_of_spl(:,1)+ mdiff(:,1); flipdim(med_of_spl(:,1)-mdiff(:,1),1)];
    fill([l_mod; flipdim(l_mod,1)], f2, [7 7 7]/8)
    hold on; plot(l , med_of_spl(:,1));

%    avg_y_mod = avg_y';
%    %figure
%    f3 = [avg_y_mod(:,1)+ mdiff(:,1); flipdim(avg_y_mod(:,1)-mdiff(:,1),1)];
%    fill([l_mod; flipdim(l_mod,1)], f3, [7 7 7]/8)
%    hold on; plot(l_mod , avg_y_mod(:,1));
   
% Apply GP to Top Half of Diff between Max and Min
pos_disc = med_of_spl(:,1)+ mdiff(:,1);

figure
plot(xs_test, pos_disc)

meanfunc = {'meanZero'};
covfunc = {'covSEiso'}; ell = 0.45; sf=0.3;
likfunc = {'likGauss'}; sn = 0.04;      
hyp = struct('mean', [], 'cov', log([ell; sf]), 'lik', log(sn));

hyp2 = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, xs_test, pos_disc);

for i = 1:length(xs_test)
    %xs_hood = linspace(xs_test(i) - 0.01, xs_test(i) + 0.01 , 10)';
    xs = linspace(0, B{1,5}(end,1), 100)';
    [pos_disc_mu pos_disc_var] = gp(hyp2, @infGaussLik, meanfunc, covfunc, likfunc, xs_test, pos_disc, xs);
    
%     pos_disc_mu(i) = mean(pos_disc_mu);
%     pos_disc_var(i) = mean(pos_disc_var);
end


 f = [pos_disc_mu(:,1)+2*sqrt(pos_disc_var(:,1)); flipdim(pos_disc_mu(:,1)-2*sqrt(pos_disc_var(:,1)),1)];
    fill([xs_test; flipdim(xs_test,1)], f, [7 7 7]/8)
    hold on; plot(xs_test, pos_disc_mu(:,1)); %plot(B{1,i}(:,1), B{1,i}(:,2), '+')

% Add Error to sources: i.e. mmax + variance of sources




%% Fusion of 2 GPs with Rho = 0
% Will need to be changed later to make more general. Also will later
% be used as reification section. Consider using
% if-else statement for the case of 2 functions to combine else 3 or more

% if number of models = 2
% for i = 1:2 %length(sheet)
    x_test = linspace(8E-8,B{1,1}(end,1),50);
    for k = 1:length(x_test)
        
        xs = linspace(x_test(k) - 8E-7 , x_test(k) + 8E-7, 10)'; 
    
        [mu1, var1] = gp(hyp2, @infGaussLik, meanfunc, covfunc, likfunc, B{1,1}(:,1), B{1,1}(:,2), xs);
        y1 = ppval(s(1),xs); %B{1,1}(:,2);% Fun_1(xs); %convert to account for function
        var1 = mu1.^2 + var1;
        
        [mu2, var2] = gp(hyp2, @infGaussLik, meanfunc, covfunc, likfunc, B{1,2}(:,1), B{1,2}(:,2), xs);
        y2 = ppval(s(2),xs); %B{1,2}(:,2); %convert to account for function
        var2 = mu2.^2 + var2;

        rho1 = zeros(length(xs),1); %zeros(length(B{1,i}(:,1)),1) ; %sqrt(var1) ./sqrt((y1-y2).^2 + var1)
        
        rho2 = zeros(length(xs),1);
        
        rho_bar(k,1) = mean((var1./(var1+var2)).*rho1 + ...  
                            (var2./(var1+var2)).*rho2);
        rho_bar(rho_bar > 0.99) = 0.99;
        rho_bar(k,2) = 0;
       
        combined_mean(k,1) = mean(WinklerMean(rho_bar(k,1), var1, var2, y1, y2)); %indeixng error here
        combined_var(k,1) = mean(WinklerVar(rho_bar(k,1), var1, var2));
        combined_mean(k,2) = mean(WinklerMean(rho_bar(k,2), var1, var2, y1, y2));
        combined_var(k,2) = mean(WinklerVar(rho_bar(k,2), var1, var2));
    end 
      
% end

figure
plot(x_test, combined_mean(:,1))









% for i = 1:2 %length(sheet)
%     x_test = linspace(8E-8,B{1,i}(end,1),50);
%     for k = 1:length(x_test)
%         
%         xs = linspace(x_test(k) - 8E-7 , x_test(k) + 8E-7, 10)'; 
%     
%         [mu{i} s2{i}] = gp(hyp2, @infGaussLik, meanfunc, covfunc, likfunc, B{1,i}(:,1), B{1,i}(:,2), xs);
% 
%         y1{i} = mu{1,i}(:,1);%B{1,i}(:,2)% Fun_1(xs); %convert to account for function
%         var1{i} = mu{1,i}(:,1).^2 + s2{1,i}(:,1);
% 
%         %y2 = Fun_2(xs); %convert to account for function
%         %var2 = mu{1,2}(:,1).^2 + s2{1,2}(:,1);
% 
%         rho1{i} = zeros(length(xs),1) ; %zeros(length(B{1,i}(:,1)),1) ; %sqrt(var1) ./sqrt((y1-y2).^2 + var1)
%         
%         rho_bar{1,i}(k,1) = mean((var1{1,i}(:,1)./(var1{1,i}(:,1)+var1{1,2}(:,1))).*rho1{1,1}(:,1) + ...  
%                             (var1{1,1}(:,1)./(var1{1,1}(:,1)+var1{1,2}(:,1))).*rho1{1,2}(:,1));
%         rho_bar(rho_bar > 0.99) = 0.99;
%         rho_bar{1,i}(k,2) = 0;
%        
%         combined_mean(j,1) = mean(WinklerMean(rho_bar(j,1), var1{1,1}(:,1), var1{1,2}(:,1), y1{1,1}(:,1), y1{1,2}(:,1))); %indeixng error here
%         combined_var(j,1) = mean(WinklerVar(rho_bar(j,1), var1{1,1}(:,1), var1{1,2}(:,1)))
%         combined_mean(j,2) = mean(WinklerMean(rho_bar(j,2), var1{1,1}(:,1), var1{1,2}(:,1), y1{1,1}(:,1), y1{1,2}(:,1)))
%         combined_var(j,2) = mean(WinklerVar(rho_bar(j,2), var1{1,1}(:,1),var1{1,2}(:,1)))
%     end 
%       
% end
        



%% Define Functions of models


%n=realmin

%% Non Used Code
% syms y(x)
% y(x) = piecewise(


%Vestiges from above loop:
%xq = linspace(0,0.15,725)'
%y = [0; A(:,2) ;0]
%plot(A(:,1),A(:,2),'o',xq,s)


%f = @(x) sin(x)

%function f(x)
%sin (x)
%end
%Not a good curve fit. 1st attempt
%f = fit(A(:,1),A(:,2), 'power2')
%plot(f,A(:,1),A(:,2))

% Code for Dr. Arroyave's Request
% figure
% xaxis1 = linspace(0,0.085);
% xaxis5 = linspace(0,0.1);
% plot(xaxis1, ppval(s(1),xaxis1)); hold on
% plot(xaxis5, ppval(s(5),xaxis5))
% 
% 
% for i = [1,5]
%     xs = linspace(0, B{1,i}(end,1), 100)';
%     [mu{i} s2{i}] = gp(hyp2, @infGaussLik, meanfunc, covfunc, likfunc, B{1,i}(:,1), B{1,i}(:,2), xs);
% 
%     f = [mu{1,i}(:,1)+2*sqrt(s2{1,i}(:,1)); flipdim(mu{1,i}(:,1)-2*sqrt(s2{1,i}(:,1)),1)];
%     fill([xs; flipdim(xs,1)], f, [7 7 7]/8)
%     hold on; plot(xs, mu{1,i}(:,1));
% end
%end of Code for Dr. Arroyave's Request