%% clearing
clear
close all
clc

%% load geometries
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

%% define Slice Direction
slice_direction_z=[0,0,1];
slice_direction_z=slice_direction_z/norm(slice_direction_z);

%% Create slices Parameters
x_spacing=0.5; % in mm
y_spacing=0.5; % in mm
z_spacing=0.75; % in mm

%% Edge Bounding Box Gap
pixel_bounds=[5,5,5];

%% load geoms
[GCS_to_SliceAlignTransMat,slice_bin_vol,slices_bin_vol,GCS_to_Image_TransMat,Image_to_GCS_TransMat,rot_TransMat,trans_TransMat,...
        scale_TransMat,x_slices, y_slices,z_slices,total_pix_num]=getTransformForGeoms(geoms,pixel_bounds,x_spacing,...
        y_spacing,z_spacing,slice_direction_z);


export_tiff=1;
export_nrrd=1;
export_mat=1;

all_slice_tic=tic;
for count_geom=1:length(geoms)
        geom=geoms{count_geom};
        filename_root=filenames{count_geom};
        [filenames_created,slices_bin_vol]=sliceCurrentGeom(geom,GCS_to_SliceAlignTransMat,Image_to_GCS_TransMat,x_slices, y_slices,z_slices,slice_bin_vol,...
                        slices_bin_vol,filename_root,export_tiff,export_nrrd,export_mat);
end
toc(all_slice_tic)