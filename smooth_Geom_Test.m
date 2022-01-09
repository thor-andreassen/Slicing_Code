%% clearing
clear
close all
clc

%% Load original file
geom=stlread('Femur_R.stl');
geom2.faces=geom.ConnectivityList;
geom2.vertices=geom.Points;
clear geom


[geom2.faces,geom2.vertices]=reducepatch(geom2.faces,geom2.vertices,4500);

patchfig=figure()
orig_patch=patch('Faces',geom2.faces,'Vertices',geom2.vertices,'FaceColor','r','EdgeAlpha',.76,...
        'FaceAlpha',.95);
%% create edge list
[geom2.edges]=findEdgesMesh(geom2.faces);

%% get Free Edges
[free_edges,watertight] = findFreeEdges(geom2.edges);

%% find node vertex normal
vertex_normals=zeros(size(geom2.vertices,1),3);
% pt1=zeros(size(geom2.vertices,1),3);
% pt2=zeros(size(geom2.vertices,1),3);
% hold on
for count_node=1:size(geom2.vertices)
        vertex_normals(count_node,:)=findVertexNormalFromMesh(geom2.faces,geom2.vertices,count_node);
%         pt1(count_node,:)=geom2.vertices(count_node,:);
%         pt2(count_node,:)=pt1(count_node,:)+vertex_normals(count_node,:)*5;
end

% hold on
% arrow(pt1,pt2);

%% fix current node
% num_node=175;
% nodes=geom2.vertices;
% num_neighbors=2;
% 
% [new_node_loc,disp]=shiftNodeLocation(geom2.faces,geom2.vertices,num_node,num_neighbors);
% 

num_neighbors=2;
max_iters=15;
counter=1;
tol=.001;
conv=Inf;
faces=geom2.faces;
nodes=geom2.vertices;

for count_node=1:size(nodes,1)
        [near_nodes_id{count_node},ring_nodes_layered{count_node}]=findKRingNodes(faces,count_node,num_neighbors);
end

tic
while counter<=max_iters && conv>tol
       current_order=1:size(nodes,1);
       current_order = current_order(randperm(length(current_order)));
       conv=-Inf;
       count_inner=1;
       for count_node=current_order
%                [new_node_loc,disp_err]=shiftNodeLocation(faces,nodes,count_node,num_neighbors);
               [new_node_loc,disp_err]=shiftNodeLocation(faces,nodes,count_node,...
                       near_nodes_id{count_node},ring_nodes_layered{count_node});

               if disp_err>=conv
                       conv=disp_err;
               end
               nodes(count_node,:)=new_node_loc;
               count_inner=count_inner+1;
               if mod(count_inner,100)==0
                       disp(num2str(count_inner));
               end
       end
       disp(conv)
       patchfig
       if counter==1
               deform_patch=patch('Faces',faces,'Vertices',nodes,'FaceColor','b','EdgeAlpha',.7,...
        'FaceAlpha',.95);
       else
               deform_patch.Vertices=nodes;
       end
       if counter==1
                histfig=figure();
       else
               histfig;
       end
       [~,qi_orig]=getMeshQuality(geom2.faces,geom2.vertices);
       [~,qi_deform]=getMeshQuality(faces,nodes);
       subplot(2,1,1);
       hist(qi_orig);
       subplot(2,1,2);
       hist(qi_deform);
       pause(.1);
       counter=counter+1
end
toc

