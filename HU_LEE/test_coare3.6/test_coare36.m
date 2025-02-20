clear; close all; clc
restoredefaultpath
addpath(genpath('~/Documents/MATLAB/COARE3.6/'))
% addpath('~/Documents/projects/sd-ni-wp/HU_LEE/')

%%
clc
u = 25:0.1:28; % mean wind speed accounting for the ocean current vector
zu = 3.5;
test = NaN*ones(size(u(:)));
test1 = ones(size(u(:)));
A = coare36vn_zrf_et(u(:),zu,27+randn(size(test)),zu,80+randn(size(test)),zu,...
    1000+randn(size(test)),28+randn(size(test)),500+randn(size(test)),NaN,...
    NaN,NaN,NaN,NaN,test,34+randn(size(test)),NaN,NaN,10,10,10);
% plot
close
% subplot(1,2,1)
plot(u,'bo'); hold on
plot(A(:,33),'ro')