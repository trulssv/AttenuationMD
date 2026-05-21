function data = extract_data(datapath, mode)
    

    display(strcat('Loading mode ', mode, 'data from scan ', datapath, '...'))

    switch mode
        
        case 'sino'
            
            datapath_sino = datapath;
            adf_data_sino = adfReadDataset(datapath_sino);
            data = adf_data_sino.VIEW.float32_data;
    
    
    
        case 'img'
    
            datapath_img = datapath;
            adf_data_img = adfReadDataset(datapath_img);
            data = adf_data_img.IMAGE.int16_data;
    
    
        case 'recon'
    
            datapath_recon = datapath;
            adf_data_recon = adfReadDataset(datapath_recon);
            adf_data_recon.IMAGE
            data = adf_data_recon.IMAGE.int16_data;
        case 'mega_mono'
            datapath_img = datapath;
            adf_data_img = adfReadDataset(datapath_img);
            data = adf_data_img.IMAGE.mega;
            
            sz = 512;
            n_slices = length(data);

            data = cat(2, data{:});
            data = data(1:2, :);
            data = cellfun(@(cell) cell.int16_data, data, 'UniformOutput', false);
            data = cat(3, data{:});
            data = reshape(data, [sz, sz, 2, n_slices]);
            data = permute(data, [3, 4, 1, 2]);
            
              
            

    end




end



