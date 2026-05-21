function plot_input(aB, aR)
    
    PE = aB{1};
    PVC = aB{2};
    DeltaPE = aR{1};
    DeltaPVC = aR{2};

    fs = 16;

    
    t0 = tiledlayout(2, 1, 'TileSpacing', 'tight');
    
    
    PElims = [-1, 1];
    PVClims = [-0.5, 0.5];
    
    t1 = tiledlayout(t0, 1, 2);
    t1.Layout.Tile = 1;
    nexttile(t1);
    t1.Padding = 'none';
    t1.TileSpacing = 'none';
    imshow(PE, PElims)
    set(gca, 'XTick', [], 'YTick', [])
    ylabel('PE', 'FontSize', fs)
    xlabel('PE VNC', 'FontSize', fs-4)
    title('(A)', 'FontSize', fs)
    
    nexttile(t1)
    imshow(DeltaPE, PElims)
    set(gca, 'XTick', [], 'YTick', [])
    xlabel('$\Delta$ PE', 'FontSize', fs-4, 'Interpreter', 'latex')
    title('(B)', 'FontSize', fs)
    
    
    c1 = colorbar;
    c1.Layout.Tile = 'east';
    c1.Position(1) = c1.Position(1) + 0.05;
    
    t2 = tiledlayout(t0, 1, 2);
    t2.Padding = 'none';
    t2.TileSpacing = 'none';

    t2.Layout.Tile = 2;
    nexttile(t2)
    imshow(PVC, PVClims)
    set(gca, 'XTick', [], 'YTick', [])
    ylabel('PVC', 'FontSize', fs)
    xlabel('PVC VNC', 'FontSize', fs-4)
    title('(C)', 'FontSize', fs)

    
    nexttile(t2)
    imshow(DeltaPVC, PVClims)
    set(gca, 'XTick', [], 'YTick', [])
    xlabel('$\Delta$ PVC', 'FontSize', fs-4, 'Interpreter', 'latex')
    title('(D)', 'FontSize', fs)
    
   
    c2 = colorbar;
    c2.Layout.Tile = 'east';
    c2.Position(3) = 0.8 * c2.Position(3);
    set(gcf, 'Position', [0, 0, 450, 500])
    
    % Save
    
    print('./visualization/plots/input','-depsc2');
    
    
