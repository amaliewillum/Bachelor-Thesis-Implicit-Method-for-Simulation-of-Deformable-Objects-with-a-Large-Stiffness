clear all
close all
clc


dut = load('CV_L_1st_e.mat');
h = dut.CV_L_1st_e.h;
delta = dut.CV_L_1st_e.delta;
L = linspace(1,9,9);

%%
figure(20)
surf(L,h,delta,log(delta))
hold on 
view(-90,0)

set(gca,'Zscale','log')
set(gca,'Yscale','log')

%%
clear all
close all
clc

dut = load('CV_rho_1st_e.mat');
h = dut.CV_rho_1st_e.h;
delta = dut.CV_rho_1st_e.delta;
rho = dut.CV_rho_1st_e.rho;

%%
close all

figure(200)
surf(rho,h,delta,log(delta))
hold on 
view(-90,0)

set(gca,'Zscale','log')
set(gca,'Yscale','log')

