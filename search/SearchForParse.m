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
    nrow = 3;
    
    subplot(nrow,nrow,1);
    vizMP(M,'motor')
    title(strcat('original:',num2str(scoreMP_NoRel(M,lib))));
    %title(strcat('original:',num2str(scoreMP(M,lib))));
    
    if verbose, fprintf(1,'searching for parse (BEGINS)\n'); end
    Do = SearchMoves(M,lib,verbose,fast_mode);   
    Do.disp_score();
    
    % gradient search
    if verbose, fprintf(1,'Optimizing Gradient : 1\n'); end
    Do.move_opt_grad();    
    
    
    subplot(nrow,nrow,2);
    vizMP(Do.M,'motor');
    title(strcat('opt grad 1:',num2str(scoreMP_NoRel(Do.M,lib))));
    
    % The stroke sub-ids have likely changed after optimization. If we 
    % don't update them now, we will be giving an advantage to stroke
    % flipping.
    if verbose, fprintf(1,'Moving subids : 1\n'); end
    Do.move_opt_subids();
    
    %h2 = figure;
    %set(h2,'visible','off')
    %vizMP(Do.M,'motor')
    %saveas(h2,strcat(dir,'2_subids.png'));
    %close(h2)
    
    subplot(nrow,nrow,3);
    vizMP(Do.M,'motor')
    
    title(strcat('move subids: ',num2str(scoreMP_NoRel(Do.M,lib))));
    % gradient search
    if verbose, fprintf(1,'Optimizing Gradient : 2\n'); end
    Do.move_opt_grad();
    Do.disp_score();
    
    subplot(nrow,nrow,4);
    vizMP(Do.M,'motor');
 
    title(strcat('opt grad 2:',num2str(scoreMP_NoRel(Do.M,lib))));   
    %h2 = figure;
    %set(h2,'visible','off')
    %vizMP(Do.M,'motor')
    %saveas(h2,strcat(dir,'3_grad.png'));
    %close(h2)
    
    % optimize the direction, order, and relations between strokes
    Do.move_opt_dir_order_rel();
    
    subplot(nrow,nrow,5);
    vizMP(Do.M,'motor');
    
    title(strcat('order + rel:',num2str(scoreMP(Do.M,lib))));
    %h2 = figure;
    %set(h2,'visible','off')
    %vizMP(Do.M,'motor')
    %saveas(h2,strcat(dir,'4_dir_ord_relations.png'));
    %close(h2)
    
    % run a split/merge search
    Do.move_split_merge();
    
    subplot(nrow,nrow,6);
    vizMP(Do.M,'motor');
    title(strcat('after spl/merge: ', num2str(scoreMP(Do.M,lib))));
    
    %h2 = figure;
    %set(h2,'visible','off')
    %vizMP(Do.M,'motor')
    %saveas(h2,strcat(dir,'5_merge.png'));
    %close(h2)
    
    saveas(h,strcat(dir,'opti_moves.png'));
    close(h)
    
    M = Do.M; 
    if verbose, fprintf(1,'searching for parse (ENDS)\n'); end
end
