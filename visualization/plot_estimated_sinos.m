function plot_estimated_sinos(AC_est, AC_gt)


    fs = 16;

    t0 = tiledlayout(2, 1, 'TileSpacing', 'none');
    
    
    Ilims = [0, 0.2];
    Gdlims = [0, 0.02];
    % Gdlims = [];
    t1 = tiledlayout(t0, 1, 2);
    t1.Layout.Tile = 1;
    nexttile(t1);
    t1.Padding = 'none';
    t1.TileSpacing = 'none';
    imshow(AC_gt{1}, Ilims)
    set(gca, 'XTick', [], 'YTick', [])
    ylabel('Iodine', 'FontSize', fs)
    xlabel('Ground Truth', 'FontSize', fs-4)
    title('(A)', 'FontSize', fs)
    
    nexttile(t1)
    imshow(AC_est{1}, Ilims)
    set(gca, 'XTick', [], 'YTick', [])
    xlabel('Prediction', 'FontSize', fs-4)
    title('(B)', 'FontSize', fs)

    
    c1 = colorbar;    
    c1.Layout.Tile = 'east';
    
    t2 = tiledlayout(t0, 1, 2);
    t2.Padding = 'none';
    t2.TileSpacing = 'none';

    t2.Layout.Tile = 2;
    nexttile(t2)
    imshow(AC_gt{2}, Gdlims)
    set(gca, 'XTick', [], 'YTick', [])
    ylabel('Gadolinium', 'FontSize', fs)
    xlabel('Ground Truth', 'FontSize', fs-4)
    title('(C)', 'FontSize', fs)

    
    nexttile(t2)
    imshow(AC_est{2}, Gdlims)
    set(gca, 'XTick', [], 'YTick', [])
    xlabel('Prediction', 'FontSize', fs-4)
    title('(D)', 'FontSize', fs)
    
   
    c2 = colorbar;
    c2.Layout.Tile = 'east';

    
    
    set(gcf, 'Position', [0, 0, 600, 500])
    
    print('./visualization/plots/sinos_estimated','-depsc2');
    