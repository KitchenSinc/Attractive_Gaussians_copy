clear
clc
close all
Amax = 10;
Amin = 0;
A1 = linspace(Amin+.001,Amax-.001,500);
A2 = A1;

attract_mean = 5;
attract_std_dev = 1.5;
min_sig = 0.5;
norm_const = integral(@(x) normpdf(x, attract_mean, attract_std_dev), Amin, Amax);
[mu, sig] = get_desire_params(A1, Amax, min_sig);
attract_pdf = @(x) normpdf(x, attract_mean, attract_std_dev)/norm_const;
desired_attract_pdf = cell(size(mu));
accept_pdf = desired_attract_pdf;
for i = 1:length(A1)
    R = integral( @(x) normpdf(x, mu(i), sig(i)), Amax, 2*Amax);
    L = integral( @(x) normpdf(x, mu(i), sig(i)), Amin-Amax, Amin);
    
    desired_attract_pdf{i} = @(x) normpdf(x, mu(i), sig(i))/(1-(R+L));
    accept_pdf{i} = @(x) integral(desired_attract_pdf{i}, Amin, x);
end

P = zeros(length(A1),length(A1));
P_joint = P;

for i = 1:length(A1)
    
    for j = 1:length(A2)
        P(j,i) = attract_pdf(A2(j))*desired_attract_pdf{i}(A2(j))*accept_pdf{j}(A1(i));
        P_joint(j,i) = P(j,i)*attract_pdf(A1(i));
    end
    
end

[~,idx] = max(P,[],1);
most_likely_mate_rating = A2(idx);
expected_match = A2*P./sum(P,1);

figure(1)
mesh(A1,A2,P)
hold on
contour3(A1,A2,P);
xlabel('Controlled Rating');
ylabel('Matchup Rating');
title('Conditional Likelihood of Matchup');
figure(2)
contourf(A1,A2,P)
xlabel('Controlled Rating');
ylabel('Matchup Rating');
title('Conditional Likelihood of Matchup');
figure(3)
plot(A1,most_likely_mate_rating, 'b', 'LineWidth',2);
hold on
plot(A1,expected_match,'r','LineWidth',2);
plot(A1,A1,'k','LineWidth',2);
xlabel('Controlled Rating');
ylabel('Matchup Rating');
title('Most Likely Matchup for Each Rating');
legend('Most likely match', 'Expected match','Perfect match')
figure(4)
mesh(A1,A2,P_joint);
hold on
contour3(A1,A2,P_joint);
xlabel('Controlled Rating');
ylabel('Matchup Rating');
title('Joint Likelihood of Matchup');
figure(5)
contourf(A1,A2,P_joint);
figure(6)
plot(A1,A2*P_joint./sum(P_joint,1),'b','LineWidth',2);

function [mu, sig] = get_desire_params(A, Amax, min_sig)
    mu = 3*A/4 + Amax/4;
    sig = max(1.5*(mu - A),min_sig*ones(size(mu-A)));
end