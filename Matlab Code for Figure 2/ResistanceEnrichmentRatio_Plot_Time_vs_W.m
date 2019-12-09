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
% function of time and the initial resistant fraction (w) at r = 16 (Figure
% 2B)


clear all; close all; clc;
load('metasim_b.mat'); %load all the workspace variable from running RareHeterogeneity_Main.m
% output figure dir
FigFolder = './Figures';
mkdir(FigFolder)
%% Organize the metrics estimated from simuations into a structure by the parameters (w,r)used 

dose_id = 2; %drug treated condition
r = 2; % R(r) = 16

for w = 1:length(W)
    
    %Homogeneously sensitive response
        
    dr3_ctrl(w).mean = squeeze(O(r,w).DR_stats(3).mat(1,1,:));
    dr3_ctrl(w).CI = squeeze(O(r,w).DR_stats(3).mat(1,2:3,:))';
    dr3(w).mean = squeeze(O(r,w).DR_stats(3).mat(dose_id,1,:));
    dr3(w).CI = squeeze(O(r,w).DR_stats(3).mat(dose_id,2:3,:))';
    
    gr3_ctrl(w).mean = squeeze(O(r,w).GR_stats(3).mat(1,1,:));
    gr3_ctrl(w).CI = squeeze(O(r,w).GR_stats(3).mat(1,2:3,:))';
    gr3(w).mean = squeeze(O(r,w).GR_stats(3).mat(dose_id,1,:));
    gr3(w).CI = squeeze(O(r,w).GR_stats(3).mat(dose_id,2:3,:))';

    dip3_ctrl(w).mean = squeeze(O(r,w).DIP_stats(3).mat(1,1,:));
    dip3_ctrl(w).CI = squeeze(O(r,w).DIP_stats(3).mat(1,2:3,:))';
    dip3(w).mean = squeeze(O(r,w).DIP_stats(3).mat(dose_id,1,:));
    dip3(w).CI = squeeze(O(r,w).DIP_stats(3).mat(dose_id,2:3,:))';
    
    k3_ctrl(w).death.mean = O(r,w).kmat_stats(3).death(:,1,1);
    k3_ctrl(w).death.CI = O(r,w).kmat_stats(3).death(:,2:3,1);
    k3(w).death.mean = O(r,w).kmat_stats(3).death(:,1,dose_id);
    k3(w).death.CI = O(r,w).kmat_stats(3).death(:,2:3,dose_id);
    
    k3_ctrl(w).divi.mean = O(r,w).kmat_stats(3).divi(:,1,1);
    k3_ctrl(w).divi.CI = O(r,w).kmat_stats(3).divi(:,2:3,1);
    k3(w).divi.mean = O(r,w).kmat_stats(3).divi(:,1,dose_id);
    k3(w).divi.CI = O(r,w).kmat_stats(3).divi(:,2:3,dose_id);
    
    k3_ctrl(w).net.mean = O(r,w).kmat_stats(3).net(:,1,1);
    k3_ctrl(w).net.CI = O(r,w).kmat_stats(3).net(:,2:3,1);
    k3(w).net.mean = O(r,w).kmat_stats(3).net(:,1,dose_id);
    k3(w).net.CI = O(r,w).kmat_stats(3).net(:,2:3,dose_id);
    
    
    %Heterogeneous response

    dr4_ctrl(w).mean = squeeze(O(r,w).DR_stats(4).mat(1,1,:));
    dr4_ctrl(w).CI = squeeze(O(r,w).DR_stats(4).mat(1,2:3,:))';
    dr4(w).mean = squeeze(O(r,w).DR_stats(4).mat(dose_id,1,:));
    dr4(w).CI = squeeze(O(r,w).DR_stats(4).mat(dose_id,2:3,:))';
    
    gr4_ctrl(w).mean = squeeze(O(r,w).GR_stats(4).mat(1,1,:));
    gr4_ctrl(w).CI = squeeze(O(r,w).GR_stats(4).mat(1,2:3,:))';
    gr4(w).mean = squeeze(O(r,w).GR_stats(4).mat(dose_id,1,:));
    gr4(w).CI = squeeze(O(r,w).GR_stats(4).mat(dose_id,2:3,:))';
    
        
    dip4_ctrl(w).mean = squeeze(O(r,w).DIP_stats(4).mat(1,1,:));
    dip4_ctrl(w).CI = squeeze(O(r,w).DIP_stats(4).mat(1,2:3,:))';
    dip4(w).mean = squeeze(O(r,w).DIP_stats(4).mat(dose_id,1,:));
    dip4(w).CI = squeeze(O(r,w).DIP_stats(4).mat(dose_id,2:3,:))';
    
    
    k4_ctrl(w).death.mean = O(r,w).kmat_stats(4).death(:,1,1);
    k4_ctrl(w).death.CI = O(r,w).kmat_stats(4).death(:,2:3,1);
    k4(w).death.mean = O(r,w).kmat_stats(4).death(:,1,dose_id);
    k4(w).death.CI = O(r,w).kmat_stats(4).death(:,2:3,dose_id);
    
    k4_ctrl(w).divi.mean = O(r,w).kmat_stats(4).divi(:,1,1);
    k4_ctrl(w).divi.CI = O(r,w).kmat_stats(4).divi(:,2:3,1);
    k4(w).divi.mean = O(r,w).kmat_stats(4).divi(:,1,dose_id);
    k4(w).divi.CI = O(r,w).kmat_stats(4).divi(:,2:3,dose_id);
    
    k4_ctrl(w).net.mean = O(r,w).kmat_stats(4).net(:,1,1);
    k4_ctrl(w).net.CI = O(r,w).kmat_stats(4).net(:,2:3,1);
    k4(w).net.mean = O(r,w).kmat_stats(4).net(:,1,dose_id);
    k4(w).net.CI = O(r,w).kmat_stats(4).net(:,2:3,dose_id);
    
    
    % Calculate resistance enrichment ratio using different drug response metrics
    
    % resistance enrichment ratio calculated using DIP metric
    dip_ctrldelta_ratio(w).mean = (dip4_ctrl(w).mean - dip4(w).mean)./ (dip3_ctrl(w).mean - dip3(w).mean);
    % resistance enrichment ratio calculated using GR metric
    gr_ctrldelta_ratio(w).mean = (gr4_ctrl(w).mean - gr4(w).mean)./ (gr3_ctrl(w).mean - gr3(w).mean);
    % resistance enrichment ratio calculated using viability
    dr_ctrldelta_ratio(w).mean = (dr4_ctrl(w).mean - dr4(w).mean)./ (dr3_ctrl(w).mean - dr3(w).mean);
    
    % resistance enrichment ratio calculated using k_death
    k_ctrldelta_ratio(w).death.mean = (k4_ctrl(w).death.mean - k4(w).death.mean)./(k3_ctrl(w).death.mean-k3(w).death.mean);
    % resistance enrichment ratio calculated using k_stasis
    k_ctrldelta_ratio(w).stasis.mean = (k4_ctrl(w).divi.mean - k4(w).divi.mean)./(k3_ctrl(w).divi.mean-k3(w).divi.mean);
  
    
 
end

% Organize resistance enrichment ratio for each metric into a matrix.
% Column: w, Row: time point
for w = 1:length(W)
    DIP_CTRLDELTA_RATIO(:,w) = dip_ctrldelta_ratio(w).mean;
    GR_CTRLDELTA_RATIO(:,w) = gr_ctrldelta_ratio(w).mean;
    DR_CTRLDELTA_RATIO(:,w) = dr_ctrldelta_ratio(w).mean;
    K_CTRLDELTA_RATIO.death(:,w) = k_ctrldelta_ratio(w).death.mean;
    K_CTRLDELTA_RATIO.stasis(:,w) = k_ctrldelta_ratio(w).stasis.mean;
end


%% Plot

figure(1)
tbins = O(1,1).tbins;
subplot(2,3,1)
imagesc(W,ceil(tbins),log2(GR_CTRLDELTA_RATIO),[-2 2]); hold on;
title('f_{a} (GR)'); xlabel('w'); ylabel('Time (hr)');
set(gca,'FontSize',12); colorbar

subplot(2,3,2)
imagesc(W,ceil(tbins),log2(DR_CTRLDELTA_RATIO),[-2 2]); hold on;
title('f_{a} (viability)'); xlabel('w'); ylabel('Time (hr)');
set(gca,'FontSize',12); colorbar 

subplot(2,3,3)
imagesc(W,ceil(tbins),log2(DIP_CTRLDELTA_RATIO),[-2 2]); hold on;
title('f_{a} (DIP)'); xlabel('w'); ylabel('Time (hr)');
set(gca,'FontSize',12); colorbar

subplot(2,3,4)
imagesc(W,ceil(tbins),log2(K_CTRLDELTA_RATIO.death),[-2 2]); hold on;
title('k_{death}'); xlabel('w'); ylabel('Time (hr)');
set(gca,'FontSize',12); colorbar

subplot(2,3,5)
imagesc(W,ceil(tbins),log2(K_CTRLDELTA_RATIO.stasis),[-2 2]); hold on;
title('k_{stasis}'); xlabel('w'); ylabel('Time (hr)');
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
savefig(h,'Figure_2B_Resistance_Enrichment_Ratio_W_vs_T.fig')
cd ..



