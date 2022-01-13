%% clearing
clear
close all
clc

%% Input Parameters


        %% define Slice Direction
        slice_direction_z=[0,0,1];
        slice_direction_z=slice_direction_z/norm(slice_direction_z);

        %% Create slices Parameters
        x_spacing=0.5; % in mm
        y_spacing=0.5; % in mm
        z_spacing=0.75; % in mm

        %% Edge Bounding Box Gap
        pixel_bounds=[5,5,5];

        %% load stl
load_stl=1;
if load_stl==1
%         filenames={'VHFL_Left_Bone_Femur.stl',...
%                 'VHFL_Muscle_VastusMedialis.stl',...
%                 'VHFL_Muscle_VastusLateralis.stl',...
%                 'VHFL_Muscle_VastusIntermedius.stl',...
%                 'VHFL_Muscle_RectusFemoris.stl'};
        filenames={'VHFL_Left_Bone_Femur.stl'};
        for geom_count=1:length(filenames)
                filename=char(filenames{geom_count});
                geom=stlread(char(filename));
                [geoms{geom_count}.faces,geoms{geom_count}.vertices]=reducepatch(geom.ConnectivityList,geom.Points,.1);
        end
        save('geom_data.mat');
else
        load('geom_data.mat');
end
%% Rotate Scan
if norm(slice_direction_z-[1,0,0])<.001
% slice along x direction
        Ny=cross([1,0,0],slice_direction_z)/norm(cross([1,0,0],slice_direction_z));
        Nx=cross(Ny,slice_direction_z)/norm(cross(Ny,slice_direction_z));
else
% all other directions
        Ny=cross([1,0,0],slice_direction_z)/norm(cross([1,0,0],slice_direction_z));
        Nx=cross(Ny,slice_direction_z)/norm(cross(Ny,slice_direction_z));
end
slice_dir_RotMat=[Nx',Ny',slice_direction_z'];

%% determine minimimum/maximum coordinate after rotation
global_vmin_rot=[Inf,Inf,Inf];
global_vmax_rot=[-Inf,-Inf,-Inf];
for geom_count=1:length(filenames)
        current_geom=geoms{geom_count}
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

% the following are the matrices that go from the world GCS coordinate
% system to the image space of the binary stack, as well as the reverse
% from the binary image stack to the corresponding space in world
% coordinates.
GCS_to_Image_TransMat=xyswap_TransMat*scale_TransMat*trans_TransMat*rot_TransMat;
Image_to_GCS_TransMat=inv(GCS_to_Image_TransMat);



slices_label_map_vol=uint8(zeros(total_pix_num(2),total_pix_num(1),length(z_slices)));
%% Slice Geometries
use_count_color=1;
for geom_count=1:length(filenames)
        filename=char(filenames{geom_count});
        current_geom=geoms{geom_count};
        % rotate geom
        geom_rot=current_geom;
        geom_rot.vertices=[slice_dir_RotMat*current_geom.vertices']';

        % align geom
        vmin_rot=min(geom_rot.vertices);
        geom_rot.vertices=geom_rot.vertices-global_vmin_rot;
        geom_rot.vertices=geom_rot.vertices+edge_bounds;

        slices_bin_vol=logical(zeros(total_pix_num(2),total_pix_num(1),length(z_slices)));

        % MAIN SLICING LOOP
        % slices in the z direction, assuming that the 
        % model is already in a bounded box. whereby the minimum
        % x, y, z, location is > 0, and the corresponding slice thicknesses,
        % and pixel reolution has been previously determined.


        [slice_sets,slices_loops,slices_polygons,slices_bin_vol]=...
        sliceGeomAlongZ(geom_rot,x_slices, y_slices,z_slices,slice_bin_vol,slices_bin_vol);

        if use_count_color
                output_slices_8bit=uint8(slices_bin_vol*geom_count);
                slices_label_map_vol=slices_label_map_vol+output_slices_8bit;
        else
                output_slices_8bit=uint8(slices_bin_vol*255);
        end
%         %% create images
%         slice_figure=figure('WindowState','Maximized');
%         for counti=1:size(slices_bin_vol,3)
% 
%               slice_plot=imshow(slices_bin_vol(:,:,counti),'Parent',gca);
%         %       slice_figure.set('WindowState','Maximized');
%               pause(.001);
%         end


        %% save image set as stack
        fname=[filename,'.tif'];
        
        options.color=false;
        options.compress='lzw';
        options.append=false;
        options.overwrite=true;
        options.big=false;
        res=saveastiff(output_slices_8bit,fname,options);



        %% Create NRRD File
        fname=[filename,'.nrrd'];
        ijkToLpsTransform=Image_to_GCS_TransMat;

        % note alothgout the structure contains the term ijktoLPS, the transform is
        % actually the transfrom from ijk to RAS.
        nrrd_img.pixelData=output_slices_8bit;
        nrrd_img.ijkToLpsTransform=ijkToLpsTransform;
        nrrdwrite(fname, nrrd_img);


        %% Create MATLAB mat File
        fname=[filename,'.mat'];
        save(fname,'slice_sets','slices_loops','slices_polygons',...
                'polygon_inside_mat','tree_array','slices_bin_vol',...
                'output_slices_8bit','GCS_to_Image_TransMat',...
                'Image_to_GCS_TransMat','ijkToLpsTransform',...
                'rot_TransMat','trans_TransMat','scale_TransMat',...
                'xyswap_TransMat','pixel_bounds','edge_bounds',...
                'slice_direction_z','global_vmin_rot','global_vmax_rot',...
                'current_geom','filename',...
                'x_slices','y_slices','z_slices','x_spacing','y_spacing','z_spacing')
end


        %% Create NRRD File
        fname=['Labelmap_Test','.nrrd'];
        ijkToLpsTransform=Image_to_GCS_TransMat;

        % note alothgout the structure contains the term ijktoLPS, the transform is
        % actually the transfrom from ijk to RAS.
        nrrd_img.pixelData=slices_label_map_vol;
        nrrd_img.ijkToLpsTransform=ijkToLpsTransform;
        nrrdwrite(fname, nrrd_img);

        fname=['Labelmap_Test','.mat'];
        save(fname,'slices_label_map_vol');