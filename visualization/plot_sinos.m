function plot_sinos(AB, AR, rotAngles, gamma)
    
    PE_sino = AB{1};
    PVC_sino = AB{2};
    DeltaPE_sino = AR{1};
    DeltaPVC_sino = AR{2};

    fs = 16;
    clf
    subplot(2, 2, 1)
    imshow(PE_sino, [], 'Xdata', rotAngles, 'Ydata', gamma)
    h1 = gca; h1.Visible = 'On';
    colorbar
    title('(A)', 'FontSize', fs)
    axis normal

    xlabel('$\theta$ (degrees)', 'interpreter', 'latex', 'FontSize', fs)
    ylabel('$\gamma$ (degrees)', 'interpreter', 'latex', 'FontSize', fs)

    subplot(2, 2, 2)
    imshow(DeltaPE_sino, [], 'Xdata', rotAngles, 'Ydata', gamma)
    h2 = gca; h2.Visible = 'On';
    colorbar
    title('(B)', 'FontSize', fs)
    axis normal

    xlabel('$\theta$ (degrees)', 'interpreter', 'latex', 'FontSize', fs)
    ylabel('$\gamma$ (degrees)', 'interpreter', 'latex', 'FontSize', fs)
    
    subplot(2, 2, 3)
    imshow(PVC_sino, [], 'Xdata', rotAngles, 'Ydata', gamma)
    h3 = gca; h3.Visible = 'On';
    colorbar
    title('(C)', 'FontSize', fs)
    axis normal

    xlabel('$\theta$ (degrees)', 'interpreter', 'latex', 'FontSize', fs)
    ylabel('$\gamma$ (degrees)', 'interpreter', 'latex', 'FontSize', fs)
    
    subplot(2, 2, 4)
    imshow(DeltaPVC_sino, [], 'Xdata', rotAngles, 'Ydata', gamma)
    h4 = gca; h4.Visible = 'On';
    colorbar
    title('(D)', 'FontSize', fs)
    axis normal

    xlabel('$\theta$ (degrees)', 'interpreter', 'latex', 'FontSize', fs)
    ylabel('$\gamma$ (degrees)', 'interpreter', 'latex', 'FontSize', fs)
    set(gcf, 'Position', [0, 0, 600, 500])
    
    print('./visualization/plots/sinos_input','-depsc2');