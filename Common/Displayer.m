classdef Displayer
	%UNTITLED4 Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
	end
	
	methods (Static)

		function plotWeights(obj,type,sol)
			
			if nargin<3;	sol=[];		end;
			
			if strcmp(type,'solution')
				axis_T2		= obj.operator.axisT2;
				weights		= sol;
				str2plot	= 'Fitted'; 
				
			elseif strcmp(type,'expmnt')
				axis_T2		= obj.DBprop.axis_T2_full;
				weights		= obj.experimental.weights;
				str2plot	= 'Experimental';
				
			end
			
			figure();
			semilogx(axis_T2,weights);grid;
			xlabel('Time [mSec]');ylabel('[AU]');
			title([str2plot,' EMC Weights']);
			
		end
		function plotWeightsComparison(Operator,Experimental,sol)
			
			axis_T2_expmnt		= Operator.DBprop.axis_T2_full;
			weights				= Experimental.experimental.weights;
			
			axis_T2_sol			= Operator.operator.axisT2;
			weightsSol			= sol;
			
			figure();
			semilogx(axis_T2_expmnt,weights,axis_T2_sol,weightsSol);grid;
			xlabel('Time [mSec]');ylabel('[AU]');
			title([' EMC Weights']);
			legend('Experimental','Reconstructed')
		end
		function plotEMC(obj,type,sol)
			
			if nargin<3;	sol=[];		end;
			TE_arr	= obj.DBprop.TE_arr;
			
			if strcmp(type,'expmnt')
				emc	= obj.experimental.e;
			elseif strcmp(type,'solution')
				emc = sol;
			end
			
			figure();
			plot(TE_arr,emc);grid;
			xlabel('Time [mSec]');ylabel('[AU]');
		end
		function plotEMCcomparison(obj,sol)
			
			TE_arr	= obj.DBprop.TE_arr;
			emcExpmnt	= obj.experimental.e;
			emcReconstructed = sol;
			emcDiff=emcExpmnt-emcReconstructed;
			sad=sum(abs(emcDiff)); % sum of absolute differences	
			
			figure();
			subplot(2,1,1);plot(TE_arr,emcExpmnt,'o-',TE_arr,emcReconstructed,'d-');grid; 
			xlabel('Time [mSec]');ylabel('[AU]');
			title('EMC'); legend('Experimental','Reconstructed');
			subplot(2,1,2);plot(TE_arr,emcDiff,'o-',TE_arr,0.*ones(size(emcDiff,1),1));grid;
			xlabel('Time [mSec]');ylabel('[AU]');
			title(sprintf('EMC difference (Experimental-Reconstructed) %.4f',sad));
			
		end
		function plotEMCcurves(obj,type)
			
			TE_arr = obj.DBprop.TE_arr;
			
			if strcmp(type,'forFit')
				str2plot	='Fitting Operator';
				E			= obj.operator.E;
			elseif strcmp(type,'expmnt')
				str2plot	='Experimental';
				E			= obj.experimental.E_noisy;
			end
			repmat(TE_arr,size(E,1),1);
			figure();plot((repmat(TE_arr,size(E,1),1))',E','.-');grid;
			xlabel('Time [mSec]');ylabel('[AU]');
			title(sprintf([str2plot,' EMC Curves (%.f #)'], size(E,1)));
		end
		
	end
	
end

