function [scoreF,score0] = mcmc_fit_type(Minit,lib,list_sid,verbose)
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

    ps = defaultps;
    nsamp = ps.mcmc.nsamp_type_chain;
    ncpt = lib.ncpt;

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
        fprintf(1,'no relations found')
        all_R = cache_enum_all_relations(lib,M);
        argmax_relations(lib,M,all_R,list_sid);
    end

    % add R.eval_spot_type if it is missing from relations --> needed for MCMC
    has_spots = M.has_eval_types(list_sid);
    if has_spots
        fprintf(1,'\neither no mid relations or all have eval spots\n');
    else
        fprintf(1,'\nsome mid relations dont have eval spots\n');
        for sid=list_sid
            if strcmp(M.S{sid}.R.type,'mid')
                if isempty(M.S{sid}.R.eval_spot_type)
                    [~,lb,ub]=bspline_gen_s(ncpt,1);
                    M.S{sid}.R.eval_spot_type=lb + rand * (ub-lb);
                end
            end
        end
    end

    % initial score
    score0 = scoreMP(M,lib,'strokes',list_sid,'type',true,'token',true,'image',true);

    % run MCMC chain for tokens
    for is = 1:nsamp
        mcmc_iter_token(MH,M,lib,list_sid);
        samples{is} = M.copy();
        curr_score = scoreMP(M,lib,'strokes',list_sid,'type',true,'token',true,'image',true);
        samples_score(is) = curr_score;
        assert(~isinf(samples_score(is)));
    end
    
    [maxscore,idx] = max(samples_score);
    best_M = samples{idx};
    final_score = scoreMP(best_M,lib,'strokes',list_sid,'type',true,'token',true,'image',true);
    fprintf(1,'\nmax score %d',maxscore);
    fprintf(1,'\nmax score 2 %d',final_score);
    thetaF = model_to_vec_fit_type(best_M,list_sid);
    refill(thetaF,Minit,list_sid);
    if ~has_relations
        scoreF = scoreMP_NoRel(Minit,lib,'strokes',list_sid,'type',true,'token',true,'image',true);
    else
        scoreF = scoreMP(Minit,lib,'strokes',list_sid,'type',true,'token',true,'image',true);
    end

    assert(Minit.has_relations==has_relations);

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

function mh_flip_token(M,lib,list_sid,searchPM)
    fprintf(1,'tryna flip a token');
    Q = M.copy()
    curr_score = scoreMP(Q,lib);
    fprintf(1,'  try flipping one stroke direction.\n');
    sid = randi([1 this.M.ns],1,1);

    UtilMP.flip_stroke(Q.S{sid});
    optimize_subids(searchPM,Q,sid);

    prop_score = scoreMP(Q,lib);
    accept = mh_accept(prop_score,curr_score);
    if accept
        fprintf(1,' flip accepted');
                % need to make sure this is legit
                M = Q.copy();
    else
        fprintf(1,' flip rejected');
    end
end

function mh_relation_token(M,searchPM)
    llvec = argmax_relations(searchPM.lib,M);
    ll = sum(llvec);
end

function mh_order_token(M,searchPM)
    optimize_order(searchPM,M)
end

function mh_split_merge_token()
        fprintf(1,' trying a split merge');
            % score current M
            curr_score = scoreMP(M,lib);
            % make a split or merge
            Q = M.copy();
            Q = SearchSplitMerge(Q,lib,verbose);
            % score new M
            prop_score = scoreMP(Q,lib);

            %% accept or reject
            accept = mh_accept(prop_score,curr_score);
            if accept
        fprintf(1,' split merge accepted');
                M = Q.copy();
        else
        fprintf(1,'split merge rejected');
            end
end

function optimize_order(searchPM,Q)
    % optimize the stroke order, where we implicitly maximize
    % over relations. But the function returns a Q where the relations are
    % not set

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
       scores(i) = optimize_relations(searchPM,QQ);
    end

    % pick the best order, and return Q without instantiating relations
    [~,windx] = randmax(scores);
    perm = P(windx,:);
    Q.S = Q.S(perm);
end
