function sol=solver(E,e,params)


nT2    = size(E,2);    % Number of T2s (nT2)


Sum_fitted_w_Eql_to_one_flag   = params.Sum_fitted_w_Eql_to_one_flag;    
Tikh_reg_flag                  = params.Tikh_reg_flag;
lambda_L1                      = params.lambda_L1;
fit_TYPE                       = params.fit_TYPE;
fit_nRep                       = params.fit_nRep;
lambda_Tikh                    = params.lambda_Tikh;
w0_Tikh                        = params.w0_Tikh;

switch fit_TYPE
    
    % Linear fitting
    case 1
        
        % Set upper / lower bounds for weights
        lb = zeros(nT2,1);
        ub = [];
        
        % Set Equality constraint: sum of all weights = 1
        if (Sum_fitted_w_Eql_to_one_flag)
            Aeq = ones(1,nT2);
            beq = 1;
        else
            Aeq = [];
            beq = [];
        end;
        
        
        options = optimoptions(@lsqlin,'MaxIter',10000);        % Set optimization options
        w_fitted = lsqlin(E,e,[],[],Aeq,beq,lb,ub,[],options);  % Perform the fit
       
        
        % Configure Tikhonov regularization
        %         if (Tikh_reg_flag)
        %             [E_sprs_chol,enew_sprs] = run_optimization_w_L1_reg_conf_Tik_reg(E,lambda_Tikh,w0_Tikh,e);
        %         else
        %             E_sprs_chol = E;
        %             %	E_full_chol = E_full;
        %             enew_sprs = e;
        %         end;
        %         
    
    % Quadratic fitting
    case 2

        % Set the initial guess
        X0 = [];
        
        for fit_idx = 1:fit_nRep
            
            % Define the Hessian matrix and f
            H = [transpose(E)*E, zeros(nT2,nT2);zeros(nT2,nT2), zeros(nT2,nT2)];
            f = [-(e)*E, lambda_L1*ones(1,nT2)];
            
            % Set upper / lower bounds for weights
            lb = zeros(2*nT2,1);
            ub = [];
            
            % Set inequality constraint
            A = [ eye(nT2) , -eye(nT2); ...
                -eye(nT2) , -eye(nT2)];
            b = zeros(2*nT2,1);
            
            % Set equality constraint: sum of all weights = 1
            if (Sum_fitted_w_Eql_to_one_flag)
                Aeq = [ones(1,nT2), zeros(1,nT2)];
                beq = 1;
            else
                Aeq = [];
                beq = [];
            end;
            
            % Set optimization options
			
            options = optimoptions(@quadprog,'MaxIter',10000);
			options.Display='off';
            
            % Perform the fit
            [wt_fitted, op_fval, op_exitflag, op_output] = quadprog(H,f,A,b,Aeq,beq,lb,ub,X0,options);
            fprintf('Finished optimization (exit code = %1.0f, n-iterations=%1.0f)\n',op_exitflag,op_output.iterations);
            
            % Set initial guess and experimental vector for next (2nd and last) iteration
            X0 = wt_fitted;

            
% 			Combine Tikhonov?
% 			if (Tikh_reg_flag)
% 				[E_sprs_chol,enew_sprs] = run_optimization_w_L1_reg_conf_Tik_reg(E,lambda_Tikh,w0_Tikh,e);
% 			else
% 				E_sprs_chol = E;
% 				enew_sprs   = e;
% 			end;

		end;
		
		w_fitted  = wt_fitted(1:nT2);
		
            %       t_fitted  = wt_fitted(nT2+1:end);
            % 		res_orig    = 0;
            % 		res_recon   = 0;
	
end;

sol=w_fitted;

end
