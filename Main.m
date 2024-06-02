


clc
clear
close all
Fun_name='F4';                    % number of test functions: 'F1' to 'F23'
SearchAgents=30;                  % number of Fishermen (population size)
Max_EFs=1000;                     % maximum number of evaluation times
[lowerbound,upperbound,dimension,fitness]=fun_info(Fun_name);                     % Object function
[Best_score,Best_pos,curve]=CFOA(SearchAgents,Max_EFs*dimension,lowerbound,upperbound,dimension,fitness);

display(['The best solution obtained by CFOA for ' [num2str(Fun_name)],'  is : ', num2str(Best_pos)]);
display(['The optimal value of the objective function found by CFOA is' [num2str(Fun_name)],'  is : ', num2str(Best_score)]);

figure
semilogy(curve,'LineWidth',2)
xlabel('Iteration');
ylabel('Best score obtained so far');
legend('CFOA')

