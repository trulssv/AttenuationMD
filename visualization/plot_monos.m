function plot_monos(a, aRef)
    
    PE = a{1};
    PVC = a{2};
    
    PE_ref = aRef{1};
    PVC_ref = aRef{2};

    E = 40;
    MuPE = 0.1*GetMu('polyethylene', E);
    MuPVC =0.1* GetMu('pvc', E);

    mono = PE * MuPE + PVC * MuPVC;
    mono = TransformHU(mono, E);

    mono_vnc = PE_ref * MuPE + PVC_ref * MuPVC;
    mono_vnc = TransformHU(mono_vnc, E);
    clf
    subplot(2,2, 1)
    imshow(mono, [])
    colorbar
    title('mono 70 keV')
    subplot(2, 2, 2)
    imshow(mono_vnc, [])
    colorbar
    title('mono 70 keV VNC')
    subplot(2, 2, 3)
    imshow(mono - mono_vnc, [])
    colorbar
    title('diff')
    subplot(2, 2, 4)
    diff_denoise = imgaussfilt(mono - mono_vnc, 10);
    imshow(diff_denoise, [])
    colorbar
    title('diff denoiseed')
    
    print('./visualization/plots/monos','-depsc2');
