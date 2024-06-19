%% Generate Measurements

%% definition for the stoping criterion
param_l2_ball.type = 'sigma';
param_l2_ball.sigma_ball = 2;


%%
q = 1;
for m = 1 : length(input_snr)
    for n = 1 : num_tests
        
        [y0, y, Nm, sigma_noise] = util_gen_measurements(x0, G, W, A, input_snr(m));
        
        [epsilon,epsilons] = util_gen_data_fidelity_bounds(y, Nm, param_l2_ball, sigma_noise);
        
        y0_t{q} = y0;
        y_t{q} = y;
        epsilont{q} = epsilon;
        epsilons_t{q} = epsilons;
        
        q = q + 1;
        
    end
end