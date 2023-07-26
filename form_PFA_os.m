function [ PFA ] = form_PFA_os( ALPHA,N,Rate )
%   ALPHA:门限因子
%   N:参考窗个数
%   Rate:比例点
 
k=ceil(N.*Rate);
PFA=gamma(N+1).*gamma(N-k+ALPHA+1)./gamma(N-k+1)./gamma(N+ALPHA+1);
 
end