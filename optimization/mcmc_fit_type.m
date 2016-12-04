function [scoreF,score0] = mcmc_fit_type(Minit,lib,list_sid,verbose,searchPM)
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

    % run Mode Hopping MCMC chain for tokens
    for is = 1:nsamp
        u_1 = rand;
        P_check_1 = 0.1;
        u_2 = rand;
        P_check_2 = 0.5;
        if u_1 > P_check_1
          % run one MCMC iteration
          mcmc_iter_token(MH,M,lib,list_sid);
          samples{is} = M.copy();
          curr_score = scoreMP(M,lib,'strokes',list_sid,'type',true,'token',true,'image',true);
          samples_score(is) = curr_score;
          assert(~isinf(samples_score(is)));
        else
          if u_2 > 0.5
            % try a flip
            %Q,searchPM.lib,list_sid,searchPM.verbose
            %Minit,lib,list_sid,verbose
            MH.mh_eval_flip_token(list_sid,M,lib,searchPM,verbose);
            samples{is} = M.copy();
            curr_score = scoreMP(M,lib,'strokes',list_sid,'type',true,'token',true,'image',true);
            samples_score(is) = curr_score;
            assert(~isinf(samples_score(is)));
	else
            % try a move_split_merge
            MH.mh_eval_split_merge_token(list_sid,M,lib,verbose);
            samples{is} = M.copy();
            curr_score = scoreMP(M,lib,'strokes',list_sid,'type',true,'token',true,'image',true);
            samples_score(is) = curr_score;
            assert(~isinf(samples_score(is)));
          end
        end
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
