%% clearing
clear
close all
clc

%% Input Parameters
        %% load stl
        filename='Femur_R.stl';
        geom=stlread(filename);
        [geom2.faces,geom2.vertices]=reducepatch(geom.ConnectivityList,geom.Points,.1);

        %% define Slice Direction
        slice_direction_z=[0,0,1];
        slice_direction_z=slice_direction_z/norm(slice_direction_z);

        %% Create slices Parameters
        x_spacing=0.5; % in mm
        y_spacing=0.5; % in mm
        z_spacing=0.75; % in mm

        %% Edge Bounding Box Gap
        pixel_bounds=[5,5,5];

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
geom_rot=geom2;
geom_rot.vertices=[slice_dir_RotMat*geom2.vertices']';

%% align
vmin_rot=min(geom_rot.vertices);
geom_rot.vertices=geom_rot.vertices-vmin_rot;

%% creating bounding box
edge_bounds=pixel_bounds.*[x_spacing,y_spacing,z_spacing];
geom_rot.vertices=geom_rot.vertices+edge_bounds;
min_bound_vertices=[0,0,0];
max_bound_vertices=max(geom_rot.vertices)+edge_bounds;
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
trans_TransMat(1:3,4)=-vmin_rot+edge_bounds;

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



%% MAIN SLICING LOOP
% slices in the z direction, assuming that the 
% model is already in a bounded box. whereby the minimum
% x, y, z, location is > 0, and the corresponding slice thicknesses,
% and pixel reolution has been previously determined.


[slice_sets,slices_loops,slices_polygons,polygon_inside_mat,tree_array,slices_bin_vol]=...
sliceGeomAlongZ(geom_rot,x_slices, y_slices,z_slices,slice_bin_vol,slices_bin_vol);

output_slices_8bit=uint8(slices_bin_vol*255);

%% create images
slice_figure=figure('WindowState','Maximized');
for counti=1:size(slices_bin_vol,3)
      
      slice_plot=imshow(slices_bin_vol(:,:,counti),'Parent',gca);
%       slice_figure.set('WindowState','Maximized');
      pause(.001);
end


%% save image set as stack
% fname='test.tif';
% 
% options.color=false;
% options.compress='lzw';
% options.append=false;
% options.overwrite=true;
% options.big=false;
% res=saveastiff(output_slices_8bit,fname,options);



%% Create NRRD File
fname='test.nrrd';
ijkToLpsTransform=Image_to_GCS_TransMat;

% note alothgout the structure contains the term ijktoLPS, the transform is
% actually the transfrom from ijk to RAS.
nrrd_img.pixelData=output_slices_8bit;
nrrd_img.ijkToLpsTransform=ijkToLpsTransform;
nrrdwrite(fname, nrrd_img);
