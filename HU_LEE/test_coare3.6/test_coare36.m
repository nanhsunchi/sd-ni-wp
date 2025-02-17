clear; close all; clc
restoredefaultpath
% addpath(genpath('~/Documents/MATLAB/COARE3.6/'))
addpath('~/Documents/projects/sd-ni-wp/HU_LEE/')

%%
clc
u = 30; % mean wind speed accounting for the ocean current vector
zu = 3.5;
test = NaN*ones(size(u(:)));
test1 = ones(size(u(:)));
A = coare36vn_zrf_et(u(:),zu,27*test1,zu,90*test1,zu,1000*test1,29*test1,500*test1,400*test1,...
    NaN,NaN,NaN,NaN,test,35*test1,test,test,10*test1,10*test1,10*test1);
% plot
close
% subplot(1,2,1)
plot(u,'bo'); hold on
plot(A(:,33),'ro')