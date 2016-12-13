function M = SearchForParse(M,lib,verbose,fast_mode,dir)
    % Search algorithm for finding a good parse
    %
    % Input
    %  M : MotorProgram
    %  verbose: display steps as we go?
    %  fast_mode: (true/false) no gradient search (by far the slowest part)

    if ~exist('verbose','var')
       verbose = false;
    end
    if ~exist('fast_mode','var')
       fast_mode = false; 
    end
    
    if isempty(M)
       return 
    end
    
    h = figure;
    set(h,'visible','off')
    sz = [500 500]; % figure size
    pos = get(h,'Position');
    pos(3:4) = sz;
    set(h,'Position',pos);
    nrow = 2;
    
    subplot(nrow,nrow,1);
    vizMP(M,'motor')
    title(strcat('original:',num2str(scoreMP_NoRel(M,lib))));
    %title(strcat('original:',num2str(scoreMP(M,lib))));
    
    if verbose, fprintf(1,'searching for parse (BEGINS)\n'); end
    Do = SearchMoves(M,lib,verbose,fast_mode);   
    Do.disp_score();

    %subplot(nrow,nrow,2);
    %vizMP(M,'motor')
    %title(strcat('relations:',num2str(scoreMP_NoRel(M,lib))));
    
    % gradient search
    if verbose, fprintf(1,'Performing MCMC : 1\n'); end
    Do.move_opt_grad(dir);    
    
    
    subplot(nrow,nrow,2);
    vizMP(Do.M,'motor');
    title(strcat('MCMC:',num2str(scoreMP(Do.M,lib))));
    
    saveas(h,strcat(dir,'opti_moves_alt.png'));
    close(h)
    
    M = Do.M; 
    if verbose, fprintf(1,'searching for parse (ENDS)\n'); end
end
