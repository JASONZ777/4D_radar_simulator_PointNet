function [ d1 ] = oscfar_alpha_binary( d1_scope,d2_precision,func,parameter)
%   d1: alpha_pre
%   d2: pha expection
%   d1_scope: alpha range
%   d1_precision: alpha preecision
%   func: form_PFA_os
%   parameter: 相关参数【类型自定】 1=ALPHA, 2=rate 3=d2
% 用于OS_CFAR门限因子的计算
    while 1
        d1=mean(d1_scope);
        d2=func(d1,parameter(1,1),parameter(1,2));   
        d2_difference=1/d2-1/parameter(1,3);
        if abs(d2_difference)<d2_precision || abs(d1_scope(1,1)-d1_scope(1,2))<0.001
            return;
        elseif d2_difference<0
            d1_scope(1,1)=d1;
        else
            d1_scope(1,2)=d1;
        end
    end
end
