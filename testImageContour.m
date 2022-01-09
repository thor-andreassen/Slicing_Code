%% clearing
clear
close all
clc

%% load data
load('Femur_R.stlsliced.mat')

%% plot data
slice_figure=figure('WindowState','Maximized');
for counti=1:size(slices_bin_vol,3)
      
      slice_plot=imshow(slices_bin_vol(:,:,counti),'Parent',gca);
%       slice_figure.set('WindowState','Maximized');
      pause(.001);
end

%% create contour

nodes=[];
% x_slice_flip=(x_slices-vmax_rot(1));
% y_slice_flip=(y_slices);
x_slice_flip=(x_slices);
y_slice_flip=(y_slices);
z_slice_flip=(z_slices);
for count_slice=1:size(slices_bin_vol,3)
        
        full_node_set=[]
        [B,L]=bwboundaries(slices_bin_vol(:,:,count_slice),'holes');
        hold on
        for k = 1:length(B)
           boundary = B{k};
           x_coords=x_slice_flip(boundary(:,2));
           y_coords=y_slice_flip(boundary(:,1));
           z_coords=ones(1,length(x_coords))*z_slice_flip(count_slice);
           nodes=[nodes;[x_coords',y_coords',z_coords']];
        end

end

%% create volume
nodes_shift=nodes-edge_bounds;
nodes_rotate=nodes_shift;
nodes_rotate=nodes_rotate+vmin_rot;

nodes_rotate=[(slice_dir_RotMat)\nodes_rotate']';
% nodes_trans=nodes_rotate+vmin_orig;
nodes_trans=nodes_rotate;

%% Plot scenarios
figure()
patch('Faces',geom.ConnectivityList,'Vertices',geom.Points,...
        'FaceColor','r','EdgeAlpha',.2)


% geom2.vertices=geom.Points;
% geom2.faces=geom.ConnectivityList;
% geom2.vertices=geom2.vertices-vmin_orig;
% patch('Faces',geom2.faces,'Vertices',geom2.vertices,...
%         'FaceColor','r','EdgeAlpha',.2)

hold on
plot3(nodes_trans(:,1),nodes_trans(:,2),nodes_trans(:,3),'bx');