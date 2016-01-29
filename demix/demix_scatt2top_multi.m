function [speech_list, lvl_1_list] = demix_scatt2top_multi(mix, Dnmf1_list, Dnmf2_list, stds1, stds2, epsf, filts, options, param1, param2, Npad)
    % mix should be a column vector 
    % speech_list is reconstructed sounds based on 2nd level (row vectors)
    % lvl_1_list is reconstructed sounds based on 1st level (row vectors)
    
    Dnmf1_all = [];
    Dnmf2_all = [];
    for i=1:length(Dnmf1_list)
        Dnmf1_all = [Dnmf1_all Dnmf1_list{i}];
        Dnmf2_all = [Dnmf2_all Dnmf2_list{i}];
    end
    
	%%% we will use a simpler algo: 
	%%% 1. compute |W1 x|
	[S2, S1, P, Plp, Php, scratch] = audioscatt_fwd_haar(pad_mirror(mix',Npad), filts, options);
	
    display('finished audioscatt_fwd_haar');
    
    %% This is for 2nd level scattering
    
	[S2r,norm1] = renorm_spect_data(S2, stds2, epsf);
	
	% run NMF on second level - this is implemented using SPAMS but other solver would work as well
    H2 = full(mexLasso(S2r,Dnmf2_all,param2));  % H is the activation
    
    display('finished H2 mexLasso');
    
    Srec2_list = {};
    for i=1:length(Dnmf2_list)
        Srec = Dnmf2_list{i} * H2(((i-1)*size(Dnmf2_list{i},2))+1:(i*size(Dnmf2_list{i},2)),:);  % reconstructed signal when using the 2nd level dictionaries
        Srec2_list{i} = Srec;
    end
    
    Srec2u_list = {};
    for i=1:length(Srec2_list)
        Srecu = unrenorm_spect_data(Srec2_list{i},stds2,norm1);
        Srec2u_list{i} = Srecu;
    end
    
    rec2_list = {};
    for i=1:length(Srec2u_list)
        rec2 = audioreconstruct2(Srec2u_list{i}, Plp, Php, filts, scratch);
        rec2_list{i} = rec2;
    end
    
    display('finished audioreconstruct2');
    
	eps = 1e-6;
    V2_ap = zeros(size(rec2_list{1}));
    for i=1:length(rec2_list)
        V2_ap = V2_ap + rec2_list{i}.^2;
    end
    V2_ap = V2_ap + eps;
    
    U2_list = {};
    for i=1:length(rec2_list)
        U2 = ((rec2_list{i}.^2)./(V2_ap)).*S1(:,1:size(V2_ap,2));
        U2_list{i} = U2;
    end
    
	%%% 2.
    Sr2_list = {}; % renormalized U
    for i=1:length(U2_list)
        Sr2 = renorm_spect_data(U2_list{i},stds1,epsf);
        Sr2_list{i} = Sr2;
    end
    
    
    
	% solve NMF on each speaker given the first solution
    H2_list = {};
    for i=1:length(Sr2_list)
        H2 = full(mexLasso(Sr2_list{i},Dnmf1_list{i},param1)); % these are the activations
        H2_list{i} = H2;
    end
    
    display('finished H2 mexLasso');
    
    rec1_list = {}; % this is the reconstructed audio based on 1st level scattering
    for i=1:length(H2_list)
        rec1 = Dnmf1_list{i}*H2_list{i};
        rec1_list{i} = rec1;
    end
    
    
	%%% 3. use mask to obtain Ui = |W1 hat{x_i} | 
	eps = 1e-6;
    V1_ap = zeros(size(rec1_list{1}));
    for i=1:length(rec1_list)
        V1_ap = V1_ap + rec1_list{i}.^2;
    end
    V1_ap = V1_ap + eps;
    
    
    U1_list = {};
    for i=1:length(rec1_list)
        U1 = ((rec1_list{i}.^2)./(V1_ap)).*S1(:,1:size(V1_ap,2));
        U1_list{i} = U1;
    end
    
    %%%%%%%%%%%%%%%%% Reconstructed Speech with 2nd level %%%%%%%%%%%%%%%%
    speech_list = {};
    for i=1:length(U1_list)
        speech = audioreconstruct1(U1_list{i}, options, filts{1}, P);
        speech_list{i} = speech;
    end
    
    display('finished audioreconstruct1');
    

    %% This is for 1st level scattering
    
	%%%% compute estimation only using first level scattering
	S1r = renorm_spect_data(S1,stds1, epsf);

        H1 =  full(mexLasso(S1r,Dnmf1_all,param1));
        
        display('finished H1 mexLasso');

        rec1_list = {}; % this is the reconstructed audio based on 1st level scattering
        for i=1:length(Dnmf1_list)
            rec1 = Dnmf1_list{i} * H1(((i-1)*size(Dnmf1_list{i},2))+1:(i*size(Dnmf1_list{i},2)),:);
            rec1_list{i} = rec1;
        end

        V1_ap = zeros(size(rec1_list{1}));
        for i=1:length(rec1_list)
            V1_ap = V1_ap + rec1_list{i}.^2;
        end
        V1_ap = V1_ap + eps;


        U1_list = {};
        for i=1:length(rec1_list)
            U1 = ((rec1_list{i}.^2)./(V1_ap)).*S1(:,1:size(V1_ap,2));
            U1_list{i} = U1;
        end


        %%%%%%%%%%%%%%%%% Reconstructed Speech with only 1st level %%%%%%%%%%%%%%%%
        lvl_1_list = {};
        for i=1:length(U1_list)
            xest = audioreconstruct1(U1_list{i}, options, filts{1}, P);
            lvl_1_list{i} = xest;
        end
        
        display('finished audioreconstruct1');

    
    %%
    
    for i=1:length(speech_list)
        speech_list{i} = speech_list{i}';
    end
    
    for i=1:length(lvl_1_list)
        lvl_1_list{i} = lvl_1_list{i}';
    end
    
