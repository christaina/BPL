function [M,scoreF,score0] = mcmc_fit_type(Minit,list_sid,searchPM,dir)
% ARGMAX_FIT_TYPE... Fit a parse to an image, with type-level parameters
% only
%
% Input
%   Minit : instance of MotorProgram class
%   list_sid: (default=all) list of stroke-ids which we want to optimize
%    while the others are frozen
%   verbose: (true/false) describe algorithm progress
%
% Output:
%  scoreF : final score
%  score0 : initial score
% 

    lib = searchPM.lib;
    verbose = false;
    ps = defaultps;
    nsamp = ps.mcmc.nsamp_type_chain;
    nsamp = 300;

    if ~exist('list_sid','var')
       list_sid = 1:M.ns; 
    end
    if ~exist('verbose','var')
       verbose = false;
    end

    MH = MCMC(false);
    samples = cell(nsamp,1);
    samples_score = nan(nsamp,1);

    M = Minit.copy();
    has_relations = M.has_relations(list_sid);

    % add relations if they are missing
    if ~has_relations
        fprintf(1,'\nno relations found')
        all_R = cache_enum_all_relations(lib,M);
        argmax_relations(lib,M,all_R,list_sid);
    end

    % add R.eval_spot_type if it is missing from relations --> needed for MCMC
    has_spots = M.has_eval_types(list_sid);
    if ~(has_spots)
        add_eval_spots(M,list_sid,lib);
    end
    now_has = M.has_eval_types(list_sid);

    % initial score
    score0 = scoreMP(M,lib,'strokes',list_sid,'type',true,'token',true,'image',true);

    % run MCMC chain for tokens

    for is = 1:nsamp
        u_1 = rand;
        P_check_1 = 0.1;
        u_2 = rand;
        now_has = M.has_eval_types(list_sid);
        if ~(now_has)
            %fprintf(1,'\n now missing spots');
            add_eval_spots(M,list_sid,lib)
        end
        if u_1 > P_check_1
            check_shapestype(M,list_sid);
            mcmc_iter_token(MH,M,lib,list_sid);
            M.clear_shapes_type();
        else 
            if u_2 > 0.9
                move_relations(M,searchPM);
            elseif (u_2 > 0.6)
                M = move_flip_all_token(M,lib,list_sid,searchPM);
            elseif (u_2 > 0.3)
                %move_order_token(M,searchPM);
                M = optimize_order(searchPM,M);
                
            else
                Q = move_split_merge_token(M,searchPM);
                fprintf(1,strcat('\nnew num strokes moved:',num2str(Q.ns)));
                M = Q.copy();
                fprintf(1,strcat('\nnew num strokes moved copy:',num2str(M.ns)));
            end
        end
        list_sid = 1:M.ns;
        curr_score = scoreMP(M,lib,'strokes',list_sid,'type',true,'token',true,'image',true);
        fprintf(1,strcat('\ni:',num2str(is),', current ll: ',num2str(curr_score)))
        samples_score(is) = curr_score;
        assert(~isinf(samples_score(is)));
        samples{is} = M.copy();
    end
    % only choose after burn-in? 
    burn_in = 100;
    [maxscore,idx] = max(samples_score(burn_in:nsamp));
    fprintf(1,strcat('\nbest sample found at idx  ',num2str(idx)));
    best_M = samples{idx+burn_in-1};
    final_score = scoreMP(best_M,lib,'strokes',1:best_M.ns,'type',true,'token',true,'image',true);
    assert(final_score==maxscore);
    fprintf(1,'\nmax score %d',maxscore);
    best_M.clear_shapes_type();
    %thetaF = model_to_vec_fit_type(best_M,1:best_M.ns);
    %refill(thetaF,Minit,list_sid);
    %Minit.S = best_M.S;
    M = best_M;
    %for is=best_M.ns;
    %    Minit.S{is}.R = best_M.S{is}.R; 
    %    Minit.S{is}.shapes_type = best_M.S{is}.shapes_type;
    %end
    
    vizSamples(samples,samples_score,nsamp,dir)
    %if ~has_relations
    %    scoreF = scoreMP_NoRel(Minit,lib,'strokes',list_sid,'type',true,'token',true,'image',true);
    %else
    scoreF = scoreMP(M,lib,'strokes',1:M.ns,'type',true,'token',true,'image',true);
    %end

    %assert(Minit.has_relations==has_relations);

end

function update_M(M_old,M_new,list_sid)
    thetaF = model_to_vec_fit_type(M_new,list_sid);
    refill(thetaF,M_old,list_sid);
    M_old.clear_relations();
    for is=list_sid
        M_old.S{is}.R = M_new.S{is}.R; 
        M_old.S{is}.shapes_type = M_new.S{is}.shapes_type;
    end

end

function check_rel(M,list_sid,lib)
    has_relations = M.has_relations(list_sid);
    if (has_relations)
        fprintf(1,'\n checking for relations and found ')
    else
    fprintf(1,'\n checking for relations and not found ')
    end

end

function refill(theta,M,list_sid)
    vec_to_model_fit_type(theta,M,list_sid);
end

% fill the MotorProgram with parameters
% and then score
function minscore = myscore_HasRel(theta,M,lib,list_sid)
    Q = M.copy(); % we don't want to modify the shared MotoProgram base
    refill(theta,Q,list_sid);
    ll = scoreMP(Q,lib,'strokes',list_sid,'type',true,'token',true,'stat',false,'image',true);    
    minscore = -ll;
end

% fill the MotorProgram with parameters
% and then optimize the relations.
% Finally, return the score "minscore"
function minscore = myscore_NoRel(theta,M,lib,list_sid,all_R)
    Q = M.copy(); % we don't want to modify the shared MotoProgram base
    refill(theta,Q,list_sid);
    argmax_relations(lib,Q,all_R,list_sid);
    ll = scoreMP(Q,lib,'strokes',list_sid,'type',true,'token',true,'stat',false,'image',true); 
    minscore = -ll;
end

function base = all_binary_strings(n)
% Generate all binary sequences of length n
% base is matrix, where rows are sequences
   base = [true; false];
   for i=2:n
      [nb,ncol] = size(base);
      base = [true(nb,1) base; false(nb,1) base];
   end
end

function M = move_flip_all_token(M,lib,list_sid,searchPM)
      if searchPM.verbose, fprintf(1,'\n move:flip stroke dirs'); end
            curr_score = scoreMP(M,searchPM.lib,'strokes',1:M.ns,'type',true,'token',true,'image',false);
            Q_flip = cell(M.ns,1);
            for sid=1:M.ns
               Q_flip{sid} = M.copy();
            Q_flip{sid}.clear_shapes_type();
               UtilMP.flip_stroke(Q_flip{sid}.S{sid});
               optimize_subids(searchPM,Q_flip{sid},sid);
               %flip_direction(searchPM,Q_flip{sid},sid);
            end

            % try all combinations of flips, 
            % with the optimal stroke order for each
            %if searchPM.verbose, fprintf(1,'\n find optimal directions/orders.'); end
            MAX_NS_ALL_PERM = 6;
            if M.ns <= MAX_NS_ALL_PERM
                bin = all_binary_strings(M.ns);
            else
                bin = rand(2^MAX_NS_ALL_PERM,M.ns)>.5;
                bin = unique(bin,'rows');
            end
            nb = size(bin,1);
            scores = zeros(nb,1);
            store_Q = cell(nb,1);
            for i=1:nb
               flip = bin(i,:);
               Q = M.copy();
               for sid=1:M.ns
                  if flip(sid) % if we flipped this stroke
                    Q.S{sid} = Q_flip{sid}.S{sid}.copy();
                  end
               end
               %optimize_order(searchPM,Q);
               optimize_relations(searchPM,Q);
               scores(i) = scoreMP(Q,searchPM.lib,'strokes',1:Q.ns,'type',true,'token',true,'image',false);
               store_Q{i} = Q;
            end

            % select the best combination of direction flips/stroke order
            [~,windx] = randmax(scores);
            % windx = argmax(scores);
            QQ = M.copy();
            for i=1:QQ.ns
                QQ.S{i}=store_Q{windx}.S{i};
            end
            prop_score = scoreMP(QQ,searchPM.lib,'strokes',1:QQ.ns,'type',true,'token',true,'image',false);

            %% accept or reject
            accept = mh_accept(prop_score,curr_score);
            if accept
                for i=1:M.ns
                    M.S{i}=store_Q{windx}.S{i};
                end
            else
                fprintf(1,'\norder move rejected');
            end
            %M = store_Q{windx}.copy();
end

function optimize_subids(searchPM,M,list_sid)
    % for all strokes in list_sid (default is all of them),
    % apply one iteration of coordinate ascent on the
    % sub-stroke ids
    if ~exist('list_sid','var')
       list_sid = 1:M.ns;
    end
    Q = M.copy();
    fprintf(1,' (optimizing SS id) ');
    sid = randi([1 M.ns],1,1);
    for sid=list_sid % each stroke
       %if searchPM.verbose, fprintf(1,'\nchoose subid for stroke %d ',sid); end
       optimize_this_subid(Q,sid,searchPM.lib,searchPM.verbose);
    end
    for i=1:M.ns
        M.S{i} = Q.S{i};
    end
    %M = Q.copy();
end

function move_relations(M,searchPM)
    Q = M.copy();
    curr_score = scoreMP(Q,searchPM.lib,'strokes',1:Q.ns,'type',true,'token',true,'image',false);
    % optimizing some relations
    llvec = argmax_relations(searchPM.lib,Q);
    ll = sum(llvec);
    fprintf(1,'move: optimize relations');
    prop_score = scoreMP(Q,searchPM.lib,'strokes',1:Q.ns,'type',true,'token',true,'image',false);
    accept = mh_accept(prop_score,curr_score);
    if accept
        %M = Q.copy();
                update_M(M,Q,1:Q.ns);
    else
        fprintf(1,'\nrelations move rejected');
    end
end

function move_order_token(M,searchPM)
    % also handles subids
    fprintf(1,'\nmove: optimize order');
    Q = M.copy();
    curr_score = scoreMP(Q,searchPM.lib,'strokes',1:Q.ns,'type',true,'token',true,'image',false);
    optimize_order(searchPM,M);
    %optimize_relations(searchPM,Q);
    %optimize_subids(searchPM,Q);
    %fprintf(1,strcat(' current aft opt ll: ',num2str(curr_score)));
    %prop_score = scoreMP(Q,searchPM.lib,'strokes',1:Q.ns,'type',true,'token',true,'image',false);

    %% accept or reject
    %accept = mh_accept(prop_score,curr_score);
    %if accept
        %M = Q.copy();
    %            update_M(M,Q,1:Q.ns);
    %else
    %    fprintf(1,'\norder move rejected');
    %end
end

function M = move_split_merge_token(M,searchPM)
        fprintf(1,' \nmove: split merge');
            % score current M
            curr_score = scoreMP(M,searchPM.lib,'strokes',1:M.ns,'type',true,'token',true,'image',false);
            % make a split or merge
            Q = M.copy();
            fprintf(1,strcat('\norig num strokes:',num2str(Q.ns)));
            % handles optimizing subids
            Q = SearchSplitMerge(Q,searchPM.lib,searchPM.verbose);
            fprintf(1,strcat('\nnew num strokes:',num2str(Q.ns)));
            % score new M
            prop_score = scoreMP(Q,searchPM.lib,'strokes',1:Q.ns,'type',true,'token',true,'image',false);

            %% accept or reject
            accept = mh_accept(prop_score,curr_score);
        if accept
                %M = Q.copy();
                %update_M(M,Q,1:Q.ns);
                M = Q;
        else
        fprintf(1,'\nsplit merge rejected');
        end
end

function M = optimize_order(searchPM,M)
    Q = M.copy();
    % optimize the stroke order, where we implicitly maximize
    % over relations. But the function returns a Q where the relations are
    % not set
    curr_score = scoreMP(Q,searchPM.lib,'strokes',1:Q.ns,'type',true,'token',true,'image',false);

    if isfield(Q.S{1},'R')
       error('cannot optimize order after relations are set');
    end

    % get all the permutations to try
    if Q.ns <= searchPM.MAX_NS_ALL_PERM
        % try all possible combinations
        P = perms(1:Q.ns);
    else
        % try a subset of the permutations
        np = factorial(searchPM.MAX_NS_ALL_PERM); %720
        P = zeros(np,Q.ns);
        for i=1:np
            P(i,:) = randperm(Q.ns);
        end
        P = unique(P,'rows');
    end

    % score all the permutations
    n = size(P,1);
    scores = zeros(n,1);
    for i=1:n
       QQ = Q.copy();
       perm = P(i,:);
       QQ.S = QQ.S(perm);
       %scores(i) = scoreMP(M,lib,'strokes',list_sid,'type',true,'token',true,'image',false);
       scores(i) = optimize_relations(searchPM,QQ);
       optimize_subids(searchPM,QQ);
    end

    % pick the best order, and return Q without instantiating relations
    [~,windx] = randmax(scores);
    perm = P(windx,:);
    Q.S = Q.S(perm);
    optimize_relations(searchPM,Q);
    optimize_subids(searchPM,Q);
    prop_score = scoreMP(Q,searchPM.lib,'strokes',1:Q.ns,'type',true,'token',true,'image',false);
    accept = mh_accept(prop_score,curr_score);
    if accept
        M = Q.copy();
    else
        fprintf(1,'\n update rejected');
    end

end

function ll = optimize_relations(searchPM,Q)
% plug in the optimal relations between the strokes
%
% ll : best log-likelihood during relation search, which includes
%      all CPDs that directly depend on the relations
    Q.clear_relations();
    llvec = argmax_relations(searchPM.lib,Q);
    ll = sum(llvec);
end

function add_eval_spots(M,list_sid,lib)
    ncpt = lib.ncpt;
    for sid=list_sid
        if strcmp(M.S{sid}.R.type,'mid')
            if isempty(M.S{sid}.R.eval_spot_type)
                [~,lb,ub]=bspline_gen_s(ncpt,1);
                M.S{sid}.R.eval_spot_type=lb + rand * (ub-lb);
            end
        end
    end
end

function check_shapestype(M,list_sid)
    for i=list_sid
        if isempty(M.S{i}.shapes_type)
            M.S{i}.shapes_type = M.S{i}.shapes_token;
        end
    end
end

function accept = mh_accept(prop_score,curr_score,g_curr_to_prop,g_prop_to_curr)
    if exist('g_curr_to_prop','var')
        lr = prop_score - curr_score + g_prop_to_curr - g_curr_to_prop;
    else
        lr = prop_score - curr_score;
    end
    fprintf(1,strcat('\noriginal score: ',num2str(curr_score)));
    fprintf(1,strcat('\nprop score: ',num2str(prop_score)));
    fprintf(1,strcat('\nscore diff: ',num2str(lr)));

    if isinf(prop_score)
       accept = false;
       return
    end

    r = min(1,exp(lr));
    assert(~isinf(curr_score));
    accept = rand <= r;
end
