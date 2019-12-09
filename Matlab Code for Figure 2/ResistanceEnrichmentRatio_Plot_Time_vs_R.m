% Author:        Natacha Comandante Lou (natacom@umich.edu)
% Copyright 2019, All Rights Reserved
%
% For paper "Phenotype-Based Probabilistic
% Analysis of Heterogeneous Responses to Cancer Drugs and Their Combination
% Efficacy" by N.Comandante-Lou, M.Khaliq, D.Venkat, M.Manikkam,
% M.Fallahi-Sichani
% ------------------------------------------------------------------------
% Description: Script for calculating resistance enrichment
% ratios using different drug response metrics and plotting them as a
% function of time and the resistance level (r) at w = 0.03 (Figure 2C)


clear all; close all; clc;
load('metasim_c.mat'); %load all the workspace variable from running RareHeterogeneity_Main.m
% output figure dir
FigFolder = './Figures';
mkdir(FigFolder)
%% Organize the metrics estimated from simuations into a structure by the parameters (w,r)used
dose_id = 2;  %drug treated condition
w = 2; % W(w) = 0.03

for r = 1:length(R)
    %Homogeneously sensitive response
    
    dip3_ctrl(r).mean = squeeze(O(r,w).DIP_stats(3).mat(1,1,:));
    dip3_ctrl(r).CI = squeeze(O(r,w).DIP_stats(3).mat(1,2:3,:))';
    dip3(r).mean = squeeze(O(r,w).DIP_stats(3).mat(dose_id,1,:));
    dip3(r).CI = squeeze(O(r,w).DIP_stats(3).mat(dose_id,2:3,:))';
    
    gr3_ctrl(r).mean = squeeze(O(r,w).GR_stats(3).mat(1,1,:));
    gr3_ctrl(r).CI = squeeze(O(r,w).GR_stats(3).mat(1,2:3,:))';
    gr3(r).mean = squeeze(O(r,w).GR_stats(3).mat(dose_id,1,:));
    gr3(r).CI = squeeze(O(r,w).GR_stats(3).mat(dose_id,2:3,:))';
    
    dr3_ctrl(r).mean = squeeze(O(r,w).DR_stats(3).mat(1,1,:));
    dr3_ctrl(r).CI = squeeze(O(r,w).DR_stats(3).mat(1,2:3,:))';
    dr3(r).mean = squeeze(O(r,w).DR_stats(3).mat(dose_id,1,:));
    dr3(r).CI = squeeze(O(r,w).DR_stats(3).mat(dose_id,2:3,:))';
    
    k3_ctrl(r).death.mean = O(r,w).kmat_stats(3).death(:,1,1);
    k3_ctrl(r).death.CI = O(r,w).kmat_stats(3).death(:,2:3,1);
    k3(r).death.mean = O(r,w).kmat_stats(3).death(:,1,dose_id);
    k3(r).death.CI = O(r,w).kmat_stats(3).death(:,2:3,dose_id);
    
    k3_ctrl(r).divi.mean = O(r,w).kmat_stats(3).divi(:,1,1);
    k3_ctrl(r).divi.CI = O(r,w).kmat_stats(3).divi(:,2:3,1);
    k3(r).divi.mean = O(r,w).kmat_stats(3).divi(:,1,dose_id);
    k3(r).divi.CI = O(r,w).kmat_stats(3).divi(:,2:3,dose_id);
    
    k3_ctrl(r).net.mean = O(r,w).kmat_stats(3).net(:,1,1);
    k3_ctrl(r).net.CI = O(r,w).kmat_stats(3).net(:,2:3,1);
    k3(r).net.mean = O(r,w).kmat_stats(3).net(:,1,dose_id);
    k3(r).net.CI = O(r,w).kmat_stats(3).net(:,2:3,dose_id);
    
    %Heterogeneous response
    dip4_ctrl(r).mean = squeeze(O(r,w).DIP_stats(4).mat(1,1,:));
    dip4_ctrl(r).CI = squeeze(O(r,w).DIP_stats(4).mat(1,2:3,:))';
    dip4(r).mean = squeeze(O(r,w).DIP_stats(4).mat(dose_id,1,:));
    dip4(r).CI = squeeze(O(r,w).DIP_stats(4).mat(dose_id,2:3,:))';
    
    gr4_ctrl(r).mean = squeeze(O(r,w).GR_stats(4).mat(1,1,:));
    gr4_ctrl(r).CI = squeeze(O(r,w).GR_stats(4).mat(1,2:3,:))';
    gr4(r).mean = squeeze(O(r,w).GR_stats(4).mat(dose_id,1,:));
    gr4(r).CI = squeeze(O(r,w).GR_stats(4).mat(dose_id,2:3,:))';
    
    dr4_ctrl(r).mean = squeeze(O(r,w).DR_stats(4).mat(1,1,:));
    dr4_ctrl(r).CI = squeeze(O(r,w).DR_stats(4).mat(1,2:3,:))';
    dr4(r).mean = squeeze(O(r,w).DR_stats(4).mat(dose_id,1,:));
    dr4(r).CI = squeeze(O(r,w).DR_stats(4).mat(dose_id,2:3,:))';
    
    k4_ctrl(r).death.mean = O(r,w).kmat_stats(4).death(:,1,1);
    k4_ctrl(r).death.CI = O(r,w).kmat_stats(4).death(:,2:3,1);
    k4(r).death.mean = O(r,w).kmat_stats(4).death(:,1,dose_id);
    k4(r).death.CI = O(r,w).kmat_stats(4).death(:,2:3,dose_id);
    
    k4_ctrl(r).divi.mean = O(r,w).kmat_stats(4).divi(:,1,1);
    k4_ctrl(r).divi.CI = O(r,w).kmat_stats(4).divi(:,2:3,1);
    k4(r).divi.mean = O(r,w).kmat_stats(4).divi(:,1,dose_id);
    k4(r).divi.CI = O(r,w).kmat_stats(4).divi(:,2:3,dose_id);
    
    k4_ctrl(r).net.mean = O(r,w).kmat_stats(4).net(:,1,1);
    k4_ctrl(r).net.CI = O(r,w).kmat_stats(4).net(:,2:3,1);
    k4(r).net.mean = O(r,w).kmat_stats(4).net(:,1,dose_id);
    k4(r).net.CI = O(r,w).kmat_stats(4).net(:,2:3,dose_id);
    
    
    
    
    % Calculate resistance enrichment ratio using different drug response metrics
    % resistance enrichment ratio calculated using DIP metric
    dip_ctrldelta_ratio(r).mean = (dip4_ctrl(r).mean - dip4(r).mean)./ (dip3_ctrl(r).mean - dip3(r).mean);
    % resistance enrichment ratio calculated using GR metric
    gr_ctrldelta_ratio(r).mean = (gr4_ctrl(r).mean - gr4(r).mean)./ (gr3_ctrl(r).mean - gr3(r).mean);
    % resistance enrichment ratio calculated using viability
    dr_ctrldelta_ratio(r).mean = (dr4_ctrl(r).mean - dr4(r).mean)./ (dr3_ctrl(r).mean - dr3(r).mean);
    % resistance enrichment ratio calculated using k_death
    k_ctrldelta_ratio(r).death.mean = (k4_ctrl(r).death.mean - k4(r).death.mean)./(k3_ctrl(r).death.mean-k3(r).death.mean);
    % resistance enrichment ratio calculated using k_stasis
    k_ctrldelta_ratio(r).stasis.mean = (k4_ctrl(r).divi.mean - k4(r).divi.mean)./(k3_ctrl(r).divi.mean-k3(r).divi.mean);
    
    
    
end

% Organize resistance enrichment ratio for each metric into a matrix.
% Column: 4, Row: time point
for r = 1:length(R)
    
    DIP_CTRLDELTA_RATIO(:,r) = dip_ctrldelta_ratio(r).mean;
    GR_CTRLDELTA_RATIO(:,r) = gr_ctrldelta_ratio(r).mean;
    DR_CTRLDELTA_RATIO(:,r) = dr_ctrldelta_ratio(r).mean;
    K_CTRLDELTA_RATIO.death(:,r) = k_ctrldelta_ratio(r).death.mean;
    K_CTRLDELTA_RATIO.stasis(:,r) = k_ctrldelta_ratio(r).stasis.mean;
    
    
end

%% Plot
figure(1)
tbins = O(1,1).tbins;
subplot(2,3,1)
imagesc(log2(R),ceil(tbins),log2(GR_CTRLDELTA_RATIO),[-2 2]); hold on;
title('f_{a} (GR)'); xlabel('log2(r)'); ylabel('Time (hr)');
set(gca,'FontSize',12); colorbar

subplot(2,3,2)
imagesc(log2(R),ceil(tbins),log2(DR_CTRLDELTA_RATIO),[-2 2]); hold on;
title('f_{a} (viability)'); xlabel('log2(r)'); ylabel('Time (hr)');
set(gca,'FontSize',12); colorbar

subplot(2,3,3)
imagesc(log2(R),ceil(tbins),log2(DIP_CTRLDELTA_RATIO),[-2 2]); hold on;
title('f_{a} (DIP)'); xlabel('log2(r)'); ylabel('Time (hr)');
set(gca,'FontSize',12); colorbar

subplot(2,3,4)
imagesc(log2(R),ceil(tbins),log2(K_CTRLDELTA_RATIO.death),[-2 2]); hold on;
title('k_{death}'); xlabel('log2(r)'); ylabel('Time (hr)');
set(gca,'FontSize',12); colorbar

subplot(2,3,5)
imagesc(log2(R),ceil(tbins),log2(K_CTRLDELTA_RATIO.stasis),[-2 2]); hold on;
title('k_{death}'); xlabel('log2(r)'); ylabel('Time (hr)');
set(gca,'FontSize',12); colorbar

%Increase the resoln of redbluecmap to make it look better
cmap = redbluecmap;
RR = cmap(:,1); G = cmap(:,2); B = cmap(:,3);
xq = linspace(1,size(cmap,1),64);
R2 = interp1(RR,xq)'; G2 = interp1(G,xq)'; B2=(interp1(B,xq))';
myrbcmap = [R2,G2,B2] ;
colormap(myrbcmap)

cd(FigFolder)
h = gcf;
savefig(h,'Figure_2C_Resistance_Enrichment_Ratio_R_vs_T.fig')
cd ..

