function [GCS_to_SliceAlignTransMat,slice_bin_vol,slices_bin_vol,GCS_to_Image_TransMat,Image_to_GCS_TransMat,rot_TransMat,trans_TransMat,...
        scale_TransMat,x_slices, y_slices,z_slices,total_pix_num]=getTransformForGeoms(geoms,pixel_bounds,x_spacing,...
        y_spacing,z_spacing,slice_direction_z)
        %% Rotate Scan
        slice_direction_z=slice_direction_z/norm(slice_direction_z);
        if norm(slice_direction_z-[1,0,0])<.001
                % slice along x direction
                Nx=cross([0,1,0],slice_direction_z)/norm(cross([0,1,0],slice_direction_z));
                Ny=cross(slice_direction_z,Nx)/norm(cross(slice_direction_z,Nx));
        else
                % all other directions
                Ny=cross([1,0,0],slice_direction_z)/norm(cross([1,0,0],slice_direction_z));
                Nx=cross(Ny,slice_direction_z)/norm(cross(Ny,slice_direction_z));
        end
        slice_dir_RotMat=[Nx',Ny',slice_direction_z'];

        %% determine minimimum/maximum coordinate after rotation
        global_vmin_rot=[Inf,Inf,Inf];
        global_vmax_rot=[-Inf,-Inf,-Inf];
        for geom_count=1:length(geoms)
                current_geom=geoms{geom_count};
                geom_rot=current_geom;
                geom_rot.vertices=[slice_dir_RotMat*current_geom.vertices']';
                vmin_rot=min(geom_rot.vertices);
                vmax_rot=max(geom_rot.vertices);
                for count_dir=1:3
                        if vmin_rot(count_dir)<=global_vmin_rot(count_dir)
                                global_vmin_rot(count_dir)=vmin_rot(count_dir);
                        end
                        if vmax_rot(count_dir)>=global_vmax_rot(count_dir)
                                global_vmax_rot(count_dir)=vmax_rot(count_dir);
                        end
                end
        end

        edge_bounds=pixel_bounds.*[x_spacing,y_spacing,z_spacing];
        min_bound_vertices=[0,0,0];
        max_bound_vertices=global_vmax_rot-global_vmin_rot+2*edge_bounds;
        x_slices=min_bound_vertices(1):x_spacing:max_bound_vertices(1);
        y_slices=min_bound_vertices(2):y_spacing:max_bound_vertices(2);
        z_slices=min_bound_vertices(3):z_spacing:max_bound_vertices(3);

        volume_pix_num=ceil((max_bound_vertices-min_bound_vertices)./[x_spacing,y_spacing,z_spacing]);
        total_pix_num=volume_pix_num;
        slice_bin_vol=false(total_pix_num(2),total_pix_num(1));
        slices_bin_vol=logical(zeros(total_pix_num(2),total_pix_num(1),length(z_slices)));

        %% Coordinate System Transformation
        % step 1 of process, rotate object to have desired slice direction axis as
        % the "z direction"
        rot_TransMat=eye(4);
        rot_TransMat(1:3,1:3)=slice_dir_RotMat;

        % step 2 of process, translate rotate object to have only positive values
        % for all coordinates. Additionally, add the translational for the bottom
        % corner spacing bounding box offset
        trans_TransMat=eye(4);
        trans_TransMat(1:3,4)=-global_vmin_rot+edge_bounds;

        % step 3 of process, scale the x,y,z coordinates to result in pixel depth
        % resolution
        scale_TransMat=eye(4);
        scale_TransMat(1,1)=1/x_spacing;
        scale_TransMat(2,2)=1/y_spacing;
        scale_TransMat(3,3)=1/z_spacing;

        % step 4 of process, flip the x and y axis for the final image creation
        xyswap_TransMat=[0,1,0,0;1,0,0,0;0,0,1,0;0,0,0,1];

        GCS_to_SliceAlignTransMat=trans_TransMat*rot_TransMat;
        
        
        % the following are the matrices that go from the world GCS coordinate
        % system to the image space of the binary stack, as well as the reverse
        % from the binary image stack to the corresponding space in world
        % coordinates.
        GCS_to_Image_TransMat=xyswap_TransMat*scale_TransMat*trans_TransMat*rot_TransMat;
        Image_to_GCS_TransMat=inv(GCS_to_Image_TransMat);
end