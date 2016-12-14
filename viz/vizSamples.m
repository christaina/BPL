function vizSamples(samples,samples_score,nsamp,dir)


    h = figure;
    set(h,'visible','off')
    %sz = [500, 500]; % figure size
    sz = 1500;
    rez = 1080;
    %set(h,'units','pixel');
    %set(h,'position',[0,0,sz,sz]);
    %set(h,'papersize',[sz,sz]);
    resolution=get(0,'ScreenPixelsPerInch'); %dont need to change anything here
    figpos=getpixelposition(h);
    set(h,'paperunits','inches','papersize',figpos(3:4)/resolution,'paperposition',[0 0 figpos(3:4)/resolution]); 
    %pos = get(h,'Position');
    %pos(3:4) = sz;
    %set(h,'Position',pos);
    nrow = 4;
    nback = 100;

    s = 1;
    %for i=(nsamp-nback+1):nsamp
    for i=1:16
        samp = randsample((nsamp-200):nsamp,1);
        samp = nsamp-(i-1);
        subplot(nrow,nrow,i);
        %title(strcat('ll:',num2str(samples_score(samp))));
        vizMP(samples{samp},'motor');
        s = s+1;
        title(strcat('ll:',num2str(samples_score(samp)),';',num2str(samp)));
    end
    print(h,fullfile(dir,'samples_alt.png'),'-dpng',['-r',num2str(rez)],'-opengl')
    %saveas(h,strcat(dir,'many_samples_3.png'));

    close(h)
end
