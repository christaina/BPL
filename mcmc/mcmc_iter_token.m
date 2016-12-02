%
% Perform an iteration of MCMC
% for re-sampling all of the token-level
% variables
%
function mcmc_iter_token(MH,M,lib,list_sid)

    if ~exist('list_sid','var')
        fprintf(1,'didnt find the sid');
        list_sid = 1:M.ns; 
    end


    % shape tokens
    for sid=list_sid
        for bid=1:M.S{sid}.nsub
            MH.mh_shape_token(sid,bid,M,lib);
        end
    end       

    % scale tokens
    for sid=list_sid
        for bid=1:M.S{sid}.nsub
            MH.mh_scale_token(sid,bid,M,lib);
        end
    end

    % position tokens
    for sid=list_sid
       MH.mh_token_position(sid,M,lib); 
    end
    
    % evaluation token spots
    for sid=list_sid
       if strcmp(M.S{sid}.R.type,'mid')
          MH.mh_eval_spot_token(sid,M,lib); 
       end        
    end

end
