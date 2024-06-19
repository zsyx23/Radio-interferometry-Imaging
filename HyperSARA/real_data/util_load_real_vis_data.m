function [y, u_, v_, nW, f, time_, pos, pixel_size] = util_load_real_vis_data(visibility_file_name, param)


speed = 299792458;
k = 1;
for m = 1 : length(visibility_file_name)
    
    m
    
    load(['./real_data/data/' visibility_file_name{m}]);
    
    
    %% Split data
    % number of overall channels
    c = size(vis,2)
    % number of spectral windows
    s = length(unique(data_id));
    
    sources = unique(field_id);
    cyg = sources(end);
    
    vis = vis(field_id==cyg,:,:);
    flaging = flaging(field_id==cyg,:,:);
    uvw = uvw(field_id==cyg,:);
    sigmas = sigmas(field_id==cyg,:);
    weights = weights(field_id==cyg,:);
    weights_ch = weights_ch(field_id==cyg,:,:);
    times = times(:,field_id==cyg)';
    data_id = data_id(:,field_id==cyg)';
    
    for i = 1 : s
        
        i
        
        Data{i} = vis(data_id==i-1,:,:);
        Flaging{i} = flaging(data_id==i-1,:,:);
        UVW{i} = uvw(data_id==i-1,:);
        Sigmas{i} = sigmas(data_id==i-1,:);
        Weights{i} = weights(data_id==i-1,:);
        Weights_ch{i} = weights_ch(data_id==i-1,:,:);
        Times{i} = times(data_id==i-1,:);
        
        %         temp = Flaging{i}(:,:,1) | Flaging{i}(:,:,4);
        %         flagI = sum(temp,2);
        %         ind{i} = (flagI==0);
        
        % sig = 0.5 * sqrt(Sigmas{i}(:,1).^2 +  Sigmas{i}(:,4).^2);
        % w_ch = 1./sig;
        
        for j = floor(c/2)
            
            temp = Flaging{i}(:,j,1) | Flaging{i}(:,j,4);
            ind = (temp==0);
            
            time_{k} = double(Times{i}(ind>0));
            
            I = (Data{i}(:,j,1) + Data{i}(:,j,4)) ./ 2;
            w_ch = sqrt(4./((1./Weights_ch{i}(:,j,1))+(1./Weights_ch{i}(:,j,4))));
            nW{k} = double(w_ch(ind>0));
            y{k} = double(I(ind>0)) .* nW{k};
            f{k} = Freqs(i,j);
            
            u = UVW{i}(:,1);
            u = u(ind>0);
            
            v = - UVW{i}(:,2);
            v = v(ind>0);
   
            
            wave = speed/f{k};
            u = u ./ wave;
            v = v ./ wave;
            bmax_k(k) = max(sqrt(u.^2 + v.^2));
            
            u_{k} = u;
            v_{k} = v;
            
            pos{k} = length(y{k});
            
            k = k + 1;
        end
        
    end
    
end

%% Setting imaging params
bmax = max(bmax_k);
theta = 1 / (2*bmax);
theta_deg = theta * 180 / pi;
theat_arcsec = theta_deg * 3600

% dl = 2.5;
% pixel_size = theat_arcsec/dl

pixel_size = 0.1871; %0.2778 % pixel size in arcsec
dl = theat_arcsec / pixel_size

%% scale the uv-coverages w.r.t the desired dl
for i = 1 : length(u_)
    v_{i} = v_{i} * pi/(bmax * dl);
    u_{i} = u_{i} * pi/(bmax * dl);
end











