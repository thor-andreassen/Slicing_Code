function [filenames_created,slices_bin_vol]=sliceCurrentGeom(geom,GCS_to_SliceAlignTransMat,Image_to_GCS_TransMat,x_slices, y_slices,z_slices,slice_bin_vol,...
                slices_bin_vol,filename_root,export_tiff,export_nrrd,export_mat)
        %% Slice Geometries
        use_count_color=0;
        
        %% rotate/align geometry
        geom_rot.vertices=GCS_to_SliceAlignTransMat*[geom.vertices';ones(1,size(geom.vertices,1))];
        geom_rot.vertices=geom_rot.vertices(1:3,:)';
        geom_rot.faces=geom.faces;
        
        % MAIN SLICING LOOP
        % slices in the z direction, assuming that the
        % model is already in a bounded box. whereby the minimum
        % x, y, z, location is > 0, and the corresponding slice thicknesses,
        % and pixel reolution has been previously determined.
        [slice_sets,slices_loops,slices_polygons,slices_bin_vol]=...
                sliceGeomAlongZ(geom_rot,x_slices, y_slices,z_slices,slice_bin_vol,slices_bin_vol*0);

        if use_count_color
                output_slices_8bit=uint8(slices_bin_vol*geom_count_color);
                slices_label_map_vol=slices_label_map_vol+output_slices_8bit;
        else
                output_slices_8bit=uint8(slices_bin_vol*255);
        end
        filenames_created={};
        file_counter=1;
        %% save image set as stack
        if export_tiff==1
                fname=[filename_root,'.tif'];
                options.color=false;
                options.compress='lzw';
                options.append=false;
                options.overwrite=true;
                options.big=false;
                res=saveastiff(output_slices_8bit,fname,options);
                filenames_created{file_counter}=fname;
                file_counter=file_counter+1;
        end


        %% Create NRRD File
        if export_nrrd==1
                fname=[filename_root,'.nrrd'];
                ijkToLpsTransform=Image_to_GCS_TransMat;

                % note alothgout the structure contains the term ijktoLPS, the transform is
                % actually the transfrom from ijk to RAS.
                nrrd_img.pixelData=output_slices_8bit;
                nrrd_img.ijkToLpsTransform=ijkToLpsTransform;
                nrrdwrite(fname, nrrd_img);
                filenames_created{file_counter}=fname;
                file_counter=file_counter+1;
        end


        %% Create MATLAB mat File
        if export_mat
                fname=[filename_root,'.mat'];
                save(fname,'slice_sets','slices_loops','slices_polygons',...
                        'slices_bin_vol','output_slices_8bit',...
                        'Image_to_GCS_TransMat','ijkToLpsTransform',...
                        'x_slices','y_slices','z_slices')
                filenames_created{file_counter}=fname;
                file_counter=file_counter+1;
        end
end