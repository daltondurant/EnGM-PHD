% plots monte carlo dots
function plot_results_MC(cfig,model_MC,truth_MC,meas_MC,est_MC,Time_store)

    n_filters = size(est_MC,2);
    names  = cell(n_filters,1);
    colors = cell(n_filters,1);

    ospa_store     = zeros(3, n_filters);
    N_store        = zeros(1, n_filters);
    G_store        = zeros(1, n_filters);

    MC = 1; % overlay this specific Monte Carlo for emphasis

    % . Plot Invividual Results
    h = [];
    for i_filters = 1:n_filters
        model = model_MC{MC,i_filters};
        truth = truth_MC;
        meas  = meas_MC{MC,i_filters};
        est   = est_MC{MC,i_filters};
        name  = model.name;
        color = model.color;

        names{i_filters}  = name;
        colors{i_filters} = color;
    
        tick_freq = 20; % number of ticks on x-axis of plots   
        strTime = num2str(round(mean(model.Time(1:tick_freq:meas.K,:),2), ...
                          3, 'significant')); % average time string array
    
        X_track = extract_tracks(truth.X,truth.total_tracks);
        
        figure()
        % This creates the 'background' axes
        set(gcf,'Color','w');

        for k=1:truth.total_tracks
            Xk= squeeze(X_track(1:3,:,k));
            plot3( Xk(1,:), Xk(2,:), Xk(3,:), '-', 'Color', [0 0 0], 'LineWidth', 1);
            hold on
            scatter3( Xk(1,1), Xk(2,1), Xk(3,1), 40, [0 0 0], "o", "filled"); % start
            hold on
            scatter3( Xk(1,end), Xk(2,end), Xk(3,end), 40, [0 0 0], "^", "filled"); % stop
        end
    
        for k = 1:meas.K
            if ~isempty(meas.C{k}) % plot clutter in state-space
                C = meas.C{k};
                s = scatter3( C(1,:), C(2,:), C(3,:),8,0.5*ones(1,3),'Marker','x');
                hold on
            end  
            if ~isempty(est.X{k}) % plot estimates
                Xk_est = est.X{k};
                plot3(Xk_est(1,:), Xk_est(2,:), Xk_est(3,:), 'LineStyle','none','Marker','.','Markersize',8,'Color',color);
                hold on
            end
        end
        Xk_est = [];
        for iMonte = 1:cfig.nMonte
            fprintf('Monte Carlo for 3D: %g / %g \n', iMonte, cfig.nMonte)
            for k = 1:meas.K
                Xk_est = cat(2, Xk_est, get_comps(est_MC{iMonte,i_filters}.X{k},[1:3]));
            end
        end
        if ~isempty(Xk_est) % plot estimates
            % for transparent data points
            %scatter3(Xk_est(1,:), Xk_est(2,:), Xk_est(3,:),8,'MarkerFaceColor',color,'Marker','o','MarkerEdgeColor','none','MarkerFaceAlpha',5/cfig.nMonte);
            % no transparency 
            scatter3(Xk_est(1,:), Xk_est(2,:), Xk_est(3,:),8,'MarkerFaceColor',color,'Marker','o','MarkerEdgeColor','none','MarkerFaceAlpha',1);
            hold on
        end

        % plot box
        c={model.range_c(1,:),model.range_c(2,:),model.range_c(3,:)};
        [c{:}]=ndgrid(c{:});
        n=length(c);
        c = reshape(cat(n+1,c{:}),[],n);
        idx = [1 5 6 2 1; 3 7 8 4 3; 5 7 3 1 5; 6 8 4 2 6; 7 8 6 5 7; 1 3 4 2 1]';
        xc = c(:,1);
        yc = c(:,2);
        zc = c(:,3);
        patch(xc(idx), yc(idx), zc(idx), 'r', 'facealpha', 0.0, 'edgecolor', 'r', 'edgealpha', 1); 
        hold on
        % pretty up
        axis equal; % keeping axis equal (no squishing)
        grid off
        ax = gca;               % get the current axis
        ax.Clipping = 'off';    % turn clipping off
        axis off;
        ax.CameraPosition  = 1.0e+03 * [1.515094243912874  -1.315469217192199   1.615107800760005];
        ax.CameraTarget    = 1.0e+02 * [1.000000000000000   0.996250267206751   2.000000000000000];
        ax.CameraUpVector  = [0 0 1];
        ax.CameraViewAngle = 11.421175255785679;
        %view([1, -1 1]);
        set(gcf,'units','inches')
        pos = get(gcf,'position');
        pos = [pos(1) 0.1*pos(2) 3.45 3.45];
        set(gcf,'position',pos)

        
        % plot OSPA errors vs. time
        figure(100); 
        if i_filters == 1
            tiles = tiledlayout(3,1,'tilespacing','tight','padding','tight');
            ax1 = nexttile(tiles,1);
            ax2 = nexttile(tiles,2);
            ax3 = nexttile(tiles,3);
        end
        ospa_c= 1e2;
        ospa_p= 2;

        ospa_vals = zeros(truth.K,3,cfig.nMonte);
        for k=1:meas.K
            for iMonte = 1:cfig.nMonte
                [ospa_vals(k,1,iMonte), ospa_vals(k,2,iMonte), ospa_vals(k,3,iMonte)]= ospa_dist(get_comps(truth_MC.X{k},[1:3]),get_comps(est_MC{iMonte,i_filters}.X{k},[1:3]),ospa_c,ospa_p);
                ospa_store(:,i_filters) = ospa_store(:,i_filters) + ospa_vals(k,:,iMonte)' ./ (cfig.nMonte*meas.K);
                N_store(:,i_filters) = N_store(:,i_filters) + (est_MC{iMonte,i_filters}.N(k))./ (cfig.nMonte*meas.K);
                G_store(:,i_filters) = G_store(:,i_filters) + (est_MC{iMonte,i_filters}.G(k))./ (cfig.nMonte*meas.K);
            end
        end
        plot(ax1, 1:meas.K, mean(ospa_vals(:,1,:),3),'Color',color)
        hold(ax1,'on')
        plot(ax2, 1:meas.K, mean(ospa_vals(:,2,:),3),'Color',color)
        hold(ax2,'on')
        plot(ax3, 1:meas.K, mean(ospa_vals(:,3,:),3),'Color',color)
        hold(ax3,'on')
        h = [h, plot(ax1, nan,nan,'LineStyle','none','Marker','.','Markersize',10,'Color',color,'DisplayName',name)];
        hold on
        
        % plot cardinality vs. time
        figure(); 
        set(gcf,'units','inches')
        pos = get(gcf,'position');
        pos = [pos(1) 0.1*pos(2) 3.45/2 3.45/2];
        set(gcf,'position',pos)
        stairs(1:meas.K,truth.total_tracks*ones(length(1:meas.K),1),'k','LineWidth',2); 
        hold on
        for iMonte =1:cfig.nMonte
            Xk_estN = est_MC{iMonte,i_filters}.N;
            scatter(1:length(Xk_estN),Xk_estN,8,'MarkerFaceColor',color,'Marker','o','MarkerEdgeColor','none','MarkerFaceAlpha',1/cfig.nMonte);
        	hold on
        end
        plot(1:meas.K,est.N,'.','Markersize',8,'Color',color);
        hold on
        set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 truth.total_tracks+10]);
        yticks([0,truth.total_tracks,10]); yticklabels({'0', '$\bf 2$', '10'});
        xlabel('Time'); ylabel('Cardinality'); xticks(1:tick_freq:meas.K); xticklabels(strTime);
        
    end % end filter loop

    % . Plot OSPA vs. time cleanup
    figure(100)
    nexttile(1);
    xlim([0 meas.K]); xticklabels(strTime); xticks(1:tick_freq:meas.K);
    ylim([0 130]); yticks([0 50 100]); 
    xlabel('(a) OSPA Distance vs. Time')
    legend(h, 'location', 'northeast', 'color', 'white', 'box', 'on', 'EdgeColor', 'white', 'NumColumns', 2)

    nexttile(2);
    xlim([0 meas.K]); xticklabels(strTime); xticks(1:tick_freq:meas.K);
    ylim([0 130]); yticks([0 50 100]); 
    xlabel('(b) OSPA Localization vs. Time')

    nexttile(3);
    xlim([0 meas.K]); xticklabels(strTime); xticks(1:tick_freq:meas.K);
    ylim([0 130]); yticks([0 50 100]); 
    xlabel('(c) OSPA Cardinality vs. Time')
    box on

    set(gcf,'units','inches')
    pos = get(gcf,'position');
    pos = [pos(1) 0.1*pos(2) 3.45 1.5*3.45];
    set(gcf,'position',pos)


    % . Plot Time averages
    figure(200)
    set(gcf,'units','inches')
    pos = get(gcf,'position');
    pos = [pos(1) 0.1*pos(2) 3.45 3.45];
    set(gcf,'position',pos)
    X = categorical(names);
    X = reordercats(X,names);
    for ii = 1:n_filters
        if ii == 1
            tiles = tiledlayout(2,1,'tilespacing','tight','padding','tight');
            ax1 = nexttile(tiles,1);
            ax2 = nexttile(tiles,2);
        end

        b = bar(ax1,X(ii),G_store(ii),'FaceColor',colors{ii}, 'edgecolor', 'none');
        xtips = b(1).XEndPoints;
        ytips = b(1).YEndPoints;
        labels = string(round(b(1).YData,3,'significant'));
        text(ax1,xtips,ytips,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom','fontsize', 8, ...
            'rotation', 30, 'backgroundcolor', 'white', 'margin', 1e-3)
        hold(ax1,'on')

        b = bar(ax2,X(ii),Time_store(ii),'FaceColor',colors{ii}, 'edgecolor', 'none');
        xtips = b(1).XEndPoints;
        ytips = b(1).YEndPoints;
        labels = string(round(b(1).YData,3,'significant'));
        text(ax2,xtips,ytips,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom','fontsize', 8, ...
            'rotation', 30, 'backgroundcolor', 'white', 'margin', 1e-3)
        hold(ax2,'on')
    end

    nexttile(1);
    xticklabels([]);
    ylim([0 500]); yticks([]);
    xlabel('(a) Average number of components')

    nexttile(2);
    ylim([0 100]); yticks([]);
    ylabel('Seconds')
    xlabel('(b) Average simulation wall-clock time')
    xtickangle(15);


    % . Plot OSPA averages
    figure(300)
    set(gcf,'units','inches')
    pos = get(gcf,'position');
    pos = [pos(1) 0.1*pos(2) 3.45 1.5*3.45];
    set(gcf,'position',pos)
    X = categorical(names);
    X = reordercats(X,names);
    for ii = 1:n_filters
        if ii == 1
            tiles = tiledlayout(3,1,'tilespacing','tight','padding','tight');
            ax1 = nexttile(tiles,1);
            ax2 = nexttile(tiles,2);
            ax3 = nexttile(tiles,3);
        end
        b = bar(ax1,ii,ospa_store(1,ii),'FaceColor',colors{ii}, 'edgecolor', 'none');
        xtips = b(1).XEndPoints; ytips = b(1).YEndPoints; labels = string(round(b(1).YData,3,'significant'));
        text(ax1, xtips,ytips,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom', 'fontsize', 8, 'rotation', 30, 'backgroundcolor', 'white', 'margin', 1e-3)
        hold(ax1,'on')

        b = bar(ax2,ii,ospa_store(2,ii),'FaceColor',colors{ii}, 'edgecolor', 'none');
        xtips = b(1).XEndPoints; ytips = b(1).YEndPoints; labels = string(round(b(1).YData,3,'significant'));
        text(ax2, xtips,ytips,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom', 'fontsize', 8, 'rotation', 30, 'backgroundcolor', 'white', 'margin', 1e-3)
        hold(ax2,'on')

        b = bar(ax3,ii,ospa_store(3,ii),'FaceColor',colors{ii}, 'edgecolor', 'none');
        xtips = b(1).XEndPoints; ytips = b(1).YEndPoints; labels = string(round(b(1).YData,3,'significant'));
        text(ax3, xtips,ytips,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom', 'fontsize', 8, 'rotation', 30, 'backgroundcolor', 'white', 'margin', 1e-3)
        hold(ax3,'on')
    end

    nexttile(1);
    xticklabels([]);
    ylim([0 110]); yticks([]);
    xlabel('(a) Time averaged OSPA distance')

    ax2 = nexttile(2);
    xticklabels([]);
    ylim([0 110]); yticks([]);
    xlabel('(b) Time averaged OSPA localization distance')

    ax3 = nexttile(3);
    ylim([0 110]); yticks([]);
    xticks([1:n_filters]); xticklabels(categorical(names)); xtickangle(15);
    xlabel('(c) Time averaged OSPA cardinality distance')
    box on

end


function X_track = extract_tracks(X,total_tracks)
    K= size(X,1); 
    x_dim= size(X{K},1); 
    k=K-1; 
    while x_dim==0 
        x_dim= size(X{k},1); 
        k= k-1; 
    end
    X_track= zeros(x_dim,K,total_tracks);

    for k=1:K
        if ~isempty(X{k})
            X_track(:,k,:)= X{k};
        end
    end
end

function Xc= get_comps(X,c)
    if isempty(X)
        Xc= [];
    else
        Xc= X(c,:);
    end
end