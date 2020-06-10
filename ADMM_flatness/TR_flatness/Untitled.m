clear all; clc;
load("data3_complete")
x_axis=ones(1,10)./[1/11:1/11:10/11]
for i =1:size(error_all_X,2)
plot(x_axis,error_all_X(:,i))
xlabel('initial conditional number for dictionary')
ylabel('error')
title('error for different initial conditional number')
legend("noise level is 0."+num2str((i-1)/2))
figure(1)
savefig(1,"D:\DL_code\simulation\TR_flatness_result\third_complete\X\noise_level_is_0."+num2str((i-1)/2)+".fig")
end




