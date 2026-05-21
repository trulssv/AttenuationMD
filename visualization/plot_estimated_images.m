function plot_estimated_images(aC_est, aC_gt)
    
    rhoGd = 7899;
    rhoI = 4944;
    
    Igt = aC_gt{1} * rhoI;
    Gdgt = aC_gt{2} * rhoGd;
    Ipred_img = aC_est{1} * rhoI;
    Gdpred_img = aC_est{2} * rhoGd;
    fs = 16;
    Idiff = abs(Igt - Ipred_img);
    Gddiff = abs(Gdgt - Gdpred_img);
    
    t0 = tiledlayout(2, 1, 'TileSpacing', 'none');
    
    
    Ilims = [0, 10];
    Gdlims = [0, 1];
    % Gdlims = [];
    t1 = tiledlayout(t0, 1, 3);
    t1.Layout.Tile = 1;
    nexttile(t1);
    t1.Padding = 'none';
    t1.TileSpacing = 'none';
    imshow(Igt, Ilims)
    set(gca, 'XTick', [], 'YTick', [])
    ylabel('Iodine', 'FontSize', fs)
    xlabel('Ground Truth', 'FontSize', fs-4)
    title('(A)', 'FontSize', fs)
    
    nexttile(t1)
    imshow(Ipred_img, Ilims)
    set(gca, 'XTick', [], 'YTick', [])
    xlabel('Prediction', 'FontSize', fs-4)
    title('(B)', 'FontSize', fs)
    
    nexttile(t1)
    imshow(Idiff, Ilims)
    set(gca, 'XTick', [], 'YTick', [])
    xlabel('Absolute Error', 'FontSize', fs-4)
    title('(C)', 'FontSize', fs)
    
    c1 = colorbar;
    c1.Layout.Tile = 'east';
    
    t2 = tiledlayout(t0, 1, 3);
    t2.Padding = 'none';
    t2.TileSpacing = 'none';

    t2.Layout.Tile = 2;
    nexttile(t2)
    imshow(Gdgt, Gdlims)
    set(gca, 'XTick', [], 'YTick', [])
    ylabel('Gadolinium', 'FontSize', fs)
    xlabel('Ground Truth', 'FontSize', fs-4)
    title('(D)', 'FontSize', fs)

    
    nexttile(t2)
    imshow(Gdpred_img, Gdlims)
    set(gca, 'XTick', [], 'YTick', [])
    xlabel('Prediction', 'FontSize', fs-4)
    title('(E)', 'FontSize', fs)
    
    
    nexttile(t2)
    imshow(Gddiff, Gdlims)
    set(gca, 'XTick', [], 'YTick', [])
    xlabel('Absolute Error', 'FontSize', fs-4)
    title('(F)', 'FontSize', fs)
    c2 = colorbar;
    c2.Layout.Tile = 'east';

    
    
    set(gcf, 'Position', [0, 0, 600, 500])
    
    print('./visualization/plots/images_estimated','-depsc2');
end