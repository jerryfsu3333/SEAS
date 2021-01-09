%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This m-file produces Tables 1, 2(in the paper), Tables 1, 2, 3(in Appendix)
% z is a 20*16 matrix, whose odd columns are the mean of the errors, and the even columns are the standard error. 
% The submitted paper only presents the mean of errors due to width limitation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Results for Model 1
clear;
tic;
z1 = zeros(20,16);
z1(1:4,:) = tab(1,100,200);
z1(5:8,:) = tab(1,100,400);
z1(9:12,:) = tab(1,100,600);
z1(13:16,:) = tab(1,200,600);
z1(17:20,:) = tab(1,400,600);
save('out_tab_M1','z1')
%
output1 = z1(:,1:2:16)
toc;
%% Results for Model 2
z2 = zeros(20,16);
z2(1:4,:) = tab(2,100,200);
z2(5:8,:) = tab(2,100,400);
z2(9:12,:) = tab(2,100,600);
z2(13:16,:) = tab(2,200,600);
z2(17:20,:) = tab(2,400,600);
save('out_tab_M2','z2')
output2 = z2(:,1:2:16)
%% Results for Model 3
z3 = zeros(20,16);
z3(1:4,:) = tab(3,100,200);
z3(5:8,:) = tab(3,100,400);
z3(9:12,:) = tab(3,100,600);
z3(13:16,:) = tab(3,200,600);
z3(17:20,:) = tab(3,400,600);
save('out_tab_M3','z3')
output3 = z3(:,1:2:16)
%% Results for Model 4
z4 = zeros(20,16);
z4(1:4,:) = tab(4,100,200);
z4(5:8,:) = tab(4,100,400);
z4(9:12,:) = tab(4,100,600);
z4(13:16,:) = tab(4,200,600);
z4(17:20,:) = tab(4,400,600);
save('out_tab_M4','z4')
output4 = z4(:,1:2:16)
%% Results for Model 5
z5 = zeros(20,16);
z5(1:4,:) = tab(5,100,200);
z5(5:8,:) = tab(5,100,400);
z5(9:12,:) = tab(5,100,600);
z5(13:16,:) = tab(5,200,600);
z5(17:20,:) = tab(5,400,600);
save('out_tab_M5','z5')
output5 = z5(:,1:2:16)
%%
toc