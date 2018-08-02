

% Descrete VS non-descrete simultation
%clc;
%clear all;
 clcl;
%clc;
% Initalizaiton
%----------------------------------------------------------------------
compName    =getComputerName;
setCD(compName);

type        = 'forFit';
sim         = Simulator;                % create a new instance of type Simulator
sim         = sim.getParams(type);		
sim         = sim.loadSimData();        
sim         = sim.simulate(type);		
axisT2      = sim.DBprop.axis_T2_full;
EMC         = sim.DBprop.emcFullNorm;
sim         = sim.getParams('exprmt');
sim         = sim.simulate('exprmt');   
axis_TE     = 1e3*sim.DBprop.TE_arr;    % [mS]


% Configuration
%----------------------------------------------------------------------
T2          = [20 45 70 140 220 310];
mixMode     = [1 2];
peakWidth   = [5 10 15 20 25 50 70 100];
nIter       = numel(T2)*numel(mixMode)*numel(peakWidth);
e           = zeros(numel(T2),numel(mixMode),numel(peakWidth),size(EMC,2));
w           = zeros(numel(T2),numel(mixMode),numel(peakWidth),size(EMC,1));

for nT2=1:numel(T2)
    sim.experimental.loc=T2(nT2);
    for nPeakW=1:numel(peakWidth)
        sim.experimental.winWidth=peakWidth(nPeakW);
        for nMode=1:numel(mixMode)
            sim.experimental.mixMode = mixMode(nMode);
            sim                      = sim.genrtWeights();	% generate vector of weights for the experimental EMC
            sim                      = sim.genrtExpEMC();	% generate experimental emc
            w(nT2,nMode,nPeakW,:)    = sim.experimental.weights;
            e(nT2,nMode,nPeakW,:)    = sim.experimental.e;
        end
    end
end


nNeighbors=5;
Neighbors_T2=(ones(nNeighbors,1)*T2)'...
            +ones(ceil(nT2),1)*(-floor(nNeighbors/2):1:floor(nNeighbors/2));

% Manually change neighbors for high T2 values
Neighbors_T2(5,[1 2 4 5])=Neighbors_T2(5,[1 2 4 5])+9;
Neighbors_T2(6,[1 2 4 5])=Neighbors_T2(6,[1 2 4 5])+20-1;
        
for i_Neighbor=1:size(Neighbors_T2,2)
    for nT2=1:numel(T2)
        sim.experimental.loc           = Neighbors_T2(nT2,i_Neighbor);
        sim.experimental.winWidth      = 0;
        sim.experimental.mixMode       = 1;
        sim                            = sim.genrtWeights();	% generate vector of weights for the experimental EMC
        sim                            = sim.genrtExpEMC();     % generate experimental emc
        w_Neighbor(nT2,i_Neighbor,:)    = sim.experimental.weights;
        e_Neighbor(nT2,i_Neighbor,:)    = sim.experimental.e;
    end
end



% Extract data for plotting
%----------------------------------------------------------------------

e_dis = squeeze(e(:,1,1,:));
for i_width=1:numel(peakWidth)
	e_gua(:,:,i_width)  = squeeze(e(:,2,i_width,:));
	er(:,:,i_width)     = 1e2*(e_gua(:,:,i_width)-e_dis)./e_dis;
% 	er(:,:,i_width)     = 1e2*(e_gua(:,:,i_width)-e_dis)./e_gua(:,:,i_width);
end

for ii=1:numel(T2)
    for jj=1:numel(peakWidth)
        avg(ii,jj)  = mean(abs(er(ii,:,jj)));
        stdd(ii,jj) = std(er(ii,:,jj));
    end
end



% 1st plot
figure;
for i_plot=1:numel(T2)
    subplot(numel(T2),1,i_plot);plot(axis_TE,abs(squeeze(er(i_plot,:,:))),'.-');grid;
    xlabel('mSec');ylabel('[%]');
    title(['EMC Difference  T_2 = ',num2str(T2(i_plot)),' mS']);
    legend(['Width ',num2str(peakWidth(1)),' mS'],['Width ',num2str(peakWidth(2)),' mS'],['Width ',num2str(peakWidth(3)),' mS'],['Width ',num2str(peakWidth(4)),' mS']);
	a=axis; axis([a(1) 240 a(3) a(4)]);
end



% 2nd plot
figure;
subplot(1,2,1);bar(avg);grid;
legend(['Width ',num2str(peakWidth(1)),' mS'],['Width ',num2str(peakWidth(2)),' mS'],['Width ',num2str(peakWidth(3)),' mS'],['Width ',num2str(peakWidth(4)),' mS']);
ylabel('[%]');
ylim([0 4]);
title('Absolute EMC Relative Difference - Average ');
barLabels={['T_2 = ',num2str(T2(1)),' mS']; ['T_2 = ',num2str(T2(2)),' mS']; ['T_2 = ',num2str(T2(3)),' mS'] };
set(gca,'xticklabel',barLabels)
subplot(1,2,2);bar(stdd);grid;
legend(['Width ',num2str(peakWidth(1)),' mS'],['Width ',num2str(peakWidth(2)),' mS'],['Width ',num2str(peakWidth(3)),' mS'],['Width ',num2str(peakWidth(4)),' mS']);
ylabel('[%]');
ylim([0 4]);
title('EMC Relative Difference - STD ');
set(gca,'xticklabel',barLabels);


% 3rd plot 
figure;
c=0;
for i_plot=1:nNeighbors
    if i_plot==ceil(nNeighbors/2)
        c=-1;
        continue;
    end
    subplot(nNeighbors-1,1,i_plot+c);
    plot(1e2*(abs((squeeze(e_Neighbor(:,i_plot,:))-e_dis)./e_dis))','.-');
    ylabel('[%]');
    legend(['T_2 = ',num2str(Neighbors_T2(i_plot,1)),' mS'],...
           ['T_2 = ',num2str(Neighbors_T2(i_plot,2)),' mS'],...
           ['T_2 = ',num2str(Neighbors_T2(i_plot,3)),' mS']);
    title(['Absolute EMC Difference of Neighbor T_2 values (Relative)']);
	grid;
end




% 4th plot
nEchos2plot=20;
x_axis2plot=ones(6,1)*(axis_TE(1:nEchos2plot));

% T2_matGau= [  squeeze(er(1,1:nEchos2plot,[1,2,3])),...
%                 squeeze(er(2,1:nEchos2plot,[2,3,4])),...
%                 squeeze(er(3,1:nEchos2plot,[2,3,4])),...
%                 squeeze(er(4,1:nEchos2plot,[2,4,6])),...
%                 squeeze(er(5,1:nEchos2plot,[4,6,8]))    ];

T2_matGau= [  squeeze(er(1,1:nEchos2plot,[1,2,3])),...
                squeeze(er(2,1:nEchos2plot,[2,3,4])),...
                squeeze(er(3,1:nEchos2plot,[3,4,5])),...
                squeeze(er(4,1:nEchos2plot,[4,5,6])),...
                squeeze(er(5,1:nEchos2plot,[6,7,8])),...
                squeeze(er(6,1:nEchos2plot,[6,7,8]))    ];


T2_matNeighbor=1e2*(abs((squeeze(e_Neighbor(:,4,1:nEchos2plot))-e_dis(:,1:nEchos2plot))./e_dis(:,1:nEchos2plot)));


close all
figure;
nGau=size(T2_matGau,2)/nT2+1;
plotVec=1:3:size(T2_matGau,2);
CM=bone(nT2);
CM(1,:)=[1,0,0];
CM(2,:)=[0,0,1];
CM(3,:)=[0,.6,0];
CM(4,:)=[0,0,0];
CM(5,:)=[.5,.1,.5];
CM(6,:)=[.1,.6,.6];

colorVec={CM(1,:);CM(2,:);CM(3,:);CM(4,:);CM(5,:);CM(6,:)};%{[0 0 1]; [1 0 0] ; [0 1 0]};

 positionVector = [0.1, 0.1, 0.45, 0.8];
    subplot('Position',positionVector)
for i_plot=1:(nT2)
    
   
    
    hNei=semilogy(axis_TE(1:nEchos2plot),abs(T2_matNeighbor(i_plot,:))','.-');
    set(hNei,'LineWidth',4)
    set(hNei,{'MarkerFaceColor'},colorVec(i_plot));
    set(hNei, {'color'}, colorVec(i_plot));
    
    hold on;
    xlabel('mSec');ylabel('[%]');
    
end




% for i_plot=1:(nT2)
%     %subplot(3,2,i_plot);
%     hold on;
%     hGau=semilogy((ones(3,1)*axis_TE(1:nEchos2plot))',abs(T2_matGau(:,plotVec(i_plot):plotVec(i_plot)+2)),'--');
%     set(hGau,{'MarkerFaceColor'},colorVec(4,1));
%     set(hGau, {'Marker'},{'o','d','*'}');
%     set(hGau, {'color'}, colorVec(i_plot));
%     set(hGau,{'MarkerFaceColor'},colorVec(i_plot));
%     set(hGau,'LineWidth',1)
%     grid on;
%     title(['T_2 = ',num2str(T2(i_plot)),' mS']);
%     %legend('Neighbor','Small Width','Midum Width','Large Width');
% end
% 
% title('Relative EMC Difference (Absolute)');
% 
% legend( ['T_2 = ',num2str(T2(1)),' mS'],...
%     ['T_2 = ',num2str(T2(2)),' mS'],...
%     ['T_2 = ',num2str(T2(3)),' mS'],...
%     ['T_2 = ',num2str(T2(4)),' mS'],...
%     ['T_2 = ',num2str(T2(5)),' mS'],...
%     ['T_2 = ',num2str(T2(6)),' mS'],...
%     'Location','southeast');
% 
% 
% width=[.6 .6 .6 .8 .8 .8];
% hight=[.1 .4 .7 .1 .4 .7 ];
% 
% 
% 
for i_plot=1:(nT2)
    
    subplot(3,2,i_plot);
   %positionVector = [width(i_plot), hight(i_plot), 0.15,0.2];
    %subplot('Position',positionVector)
    
    hNei=semilogy(axis_TE(1:nEchos2plot),abs(T2_matNeighbor(i_plot,:))','.-');
    set(hNei,'LineWidth',4)
    set(hNei,{'MarkerFaceColor'},colorVec(i_plot));
    set(hNei, {'color'}, colorVec(i_plot));
    
    hold on;
    xlabel('mSec');ylabel('[%]');
    
end

% legend( ['T_2 = ',num2str(T2(1)),' mS'],...
% 		['T_2 = ',num2str(T2(2)),' mS'],...
% 		['T_2 = ',num2str(T2(3)),' mS'],...
% 		['T_2 = ',num2str(T2(4)),' mS'],...
% 		['T_2 = ',num2str(T2(5)),' mS'],...
% 		['T_2 = ',num2str(T2(6)),' mS']);

c=1;
for i_plot=1:(nT2)
     hold on;
    subplot(3,2,i_plot);
    %positionVector = [width(i_plot), hight(i_plot), 0.15,0.2];
    %subplot('Position',positionVector)
    hGau=semilogy((ones(3,1)*axis_TE(1:nEchos2plot))',abs(T2_matGau(:,plotVec(i_plot):plotVec(i_plot)+2)),'--');
    set(hGau,{'MarkerFaceColor'},colorVec(4,1));
    set(hGau, {'Marker'},{'s','d','*'}');
    set(hGau, {'color'}, colorVec(i_plot));
    set(hGau,{'MarkerFaceColor'},colorVec(i_plot));
    set(hGau,'LineWidth',1)
    grid on;
    title(['T_2 = ',num2str(T2(i_plot)),' mS']);
	vv=peakWidth([1,2,3]+i_plot-1);
	if i_plot==(nT2) ; c=10;	end; % Neighbors of T2 values higher than 300mS are different in 10mS   
    legend(	['Neighbor T2 = ', num2str(T2(i_plot)+c) '[mS]'],...
			['Width: ',num2str(vv(1)),' [mS]'],...
			['Width: ',num2str(vv(2)),' [mS]'],...
			['Width: ',num2str(vv(3)),' [mS]'],...
		    'Location','southeast');
       
       ylim([0 1e2]);
end



