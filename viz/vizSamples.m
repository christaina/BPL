function vizSamples(samples,samples_score,nsamp,dir)


    h = figure;
    set(h,'visible','off')
    sz = [2000, 2000]; % figure size
    pos = get(h,'Position');
    pos(3:4) = sz;
    set(h,'Position',pos);
    nrow = 4;
    nback = 100;

    s = 1;
    %for i=(nsamp-nback+1):nsamp
    for i=1:16
        samp = randsample(nsamp,nsamp-100);
        subplot(nrow,nrow,s);
        s = s+1;
        vizMP(samples{samp},'motor')
        %title(strcat('ll:',num2str(samples_score(samp)),' idx:',num2str(samp)));
        title(strcat('ll:',num2str(samples_score(samp))));
    end
    saveas(h,strcat(dir,'many_samples.png'));
    close(h)
end
