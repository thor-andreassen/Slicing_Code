%% clearing
clear
close all
clc

%% load data
load('muscle_geom_orig.mat');

num_neighbors=2;
max_iters=15;
counter=1;
tol=.001;
conv=Inf;
faces=geom1.faces;
nodes=geom1.vertices;

aspects=zeros(size(faces,1),1);
skewness=zeros(size(faces,1),1);


patchfig=figure();
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
               try
                    [new_node_loc,disp_err]=shiftNodeLocation(faces,nodes,count_node,...
                           near_nodes_id{count_node},ring_nodes_layered{count_node});
               catch
                   new_node_loc=nodes(count_node,:);
               end

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
       patchfig;
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
%        [~,qi_orig]=getMeshQuality(geom1.faces,geom1.vertices);
%        [~,qi_deform]=getMeshQuality(faces,nodes);


    
        for count_face=1:size(faces,1)
            nodel=faces(count_face,:);
            face_nodes=geom1.vertices(nodel,:);
            [skewness_orig(count_face),aspects_orig(count_face)]=getMeshQuality2(face_nodes,1);

        end

         for count_face=1:size(faces,1)
            nodel=faces(count_face,:);
            face_nodes=nodes(nodel,:);
            [skewness_orig(count_face),aspects_new(count_face)]=getMeshQuality2(face_nodes,1);

        end
        
        

       clf
       cdfplot(aspects_orig);
       hold on
       cdfplot(aspects_new);
       legend('Original','Improved');
       xlim([0,5])
       pause(.1);
       counter=counter+1
end
toc