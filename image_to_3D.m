%% clearing
clear
close all
clc

%% Load transform
load('VHFL_Muscle_VastusMedialis.stl.mat','Image_to_GCS_TransMat');

%% Load original file
geom=stlread('VHFL_Muscle_VastusMedialis.stl');

%% load data
imageFullFileName='VHFL_Muscle_VastusMedialis.stl.tif';
info = imfinfo(imageFullFileName);
image_width = info(1).Width;
image_height = info(2).Height;
image_depth = length(info);

image_data=zeros(image_height,image_width,image_depth);
for count_slice = 1 : image_depth
    current_slice = imread(imageFullFileName, count_slice);
    image_data(:,:,count_slice)=current_slice;
end		


%% plot data
% slice_figure=figure('WindowState','Maximized');
% for counti=1:size(image_data,3)
%       slice_plot=imshow(image_data(:,:,counti),'Parent',gca);
%       pause(.001);
% end

%% create contour

nodes=[];
for count_slice=1:size(image_data,3)
        [B,L]=bwboundaries(image_data(:,:,count_slice),'holes');
        hold on
        for count_boundary = 1:length(B)
           boundary = B{count_boundary};
           x_coords=boundary(:,1);
           y_coords=boundary(:,2);
           z_coords=ones(1,length(x_coords))*count_slice;
           nodes=[nodes;[x_coords,y_coords,z_coords']];
        end

end

%% create volume
nodes_original=Image_to_GCS_TransMat*[nodes';ones(1,size(nodes,1))];
nodes_trans=nodes_original(1:3,:)';
nodes_trans=unique(nodes_trans,'rows');
%% Plot scenarios
figure()
patch('Faces',geom.ConnectivityList,'Vertices',geom.Points,...
        'FaceColor','r','EdgeAlpha',.2)

hold on
plot3(nodes_trans(:,1),nodes_trans(:,2),nodes_trans(:,3),'bx');

%% create Triangulation


image_surf.faces=MyCrustOpen(nodes_trans);
image_surf.vertices=nodes_trans;


[image_surf_reduce.faces,image_surf_reduce.vertices]=reducepatch(image_surf.faces,...
        image_surf.vertices,.1)
patch('Faces',image_surf_reduce.faces,'Vertices',image_surf_reduce.vertices,...
        'FaceColor','g','EdgeAlpha',.2)