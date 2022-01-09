% % %% clearing
% % clear
% % close all
% % clc
% % 
% % %% load stl
% % stl_tic=tic;
% % filename='Femur_R.stl';
% % geom=stlread(filename);
% % [geom2.faces,geom2.vertices]=reducepatch(geom.ConnectivityList,geom.Points,.1);
% % toc(stl_tic)
% % %% align
% % vmin_orig=[0,0,0];
% % vmax_orig=max(geom2.vertices);
% % 
% % % vmin_orig=min(geom2.vertices);
% % % vmax_orig=max(geom2.vertices);
% % 
% % geom2.vertices=geom2.vertices-vmin_orig;
% % vmin=min(geom2.vertices);
% % vmax=max(geom2.vertices);
% % 
% % %% define Slice Direction
% % slice_direction_z=[0,1,0];
% % slice_direction_z=slice_direction_z/norm(slice_direction_z);
% % 
% % %% Create slices Parameters
% % x_spacing=0.5; % in mm
% % y_spacing=0.5; % in mm
% % z_spacing=0.75; % in mm
% % 
% % 
% % %% Rotate Scan
% % if norm(slice_direction_z-[1,0,0])<.001
% % % slice along x direction
% %         Ny=cross([1,0,0],slice_direction_z)/norm(cross([1,0,0],slice_direction_z));
% %         Nx=cross(Ny,slice_direction_z)/norm(cross(Ny,slice_direction_z));
% % else
% % % all other directions
% %         Ny=cross([1,0,0],slice_direction_z)/norm(cross([1,0,0],slice_direction_z));
% %         Nx=cross(Ny,slice_direction_z)/norm(cross(Ny,slice_direction_z));
% % end
% % slice_dir_RotMat=[Nx',Ny',slice_direction_z'];
% % geom_rot=geom2;
% % geom_rot.vertices=[slice_dir_RotMat*geom2.vertices']';
% % 
% % %% align
% % vmin_rot=min(geom_rot.vertices);
% % vmax_rot=max(geom_rot.vertices);
% % geom_rot.vertices=geom_rot.vertices-vmin_rot;
% % 
% % %% creating bounding box
% % pixel_bounds=[5,5,5];
% % edge_bounds=pixel_bounds.*[x_spacing,y_spacing,z_spacing];
% % geom_rot.vertices=geom_rot.vertices+edge_bounds;
% % min_bound_vertices=[0,0,0];
% % max_bound_vertices=max(geom_rot.vertices)+edge_bounds;
% % x_slices=min_bound_vertices(1):x_spacing:max_bound_vertices(1);
% % y_slices=min_bound_vertices(2):y_spacing:max_bound_vertices(2);
% % z_slices=min_bound_vertices(3):z_spacing:max_bound_vertices(3);
% % 
% % volume_pix_num=ceil((max_bound_vertices-min_bound_vertices)./[x_spacing,y_spacing,z_spacing]);
% % total_pix_num=volume_pix_num;
% % slice_bin_vol=false(total_pix_num(2),total_pix_num(1));
% % slices_bin_vol=logical(zeros(total_pix_num(2),total_pix_num(1),length(z_slices)));
% % 
% % %% MAIN SLICING LOOP
% % % slices in the z direction, assuming that the 
% % % model is already in a bounded box. whereby the minimum
% % % x, y, z, location is > 0, and the corresponding slice thicknesses,
% % % and pixel reolution has been previously determined.

function [slice_sets,slices_loops,slices_polygons,polygon_inside_mat,tree_array,slices_bin_vol]=...
sliceGeomAlongZ(geom_rot,x_slices, y_slices,z_slices,slices_bin_vol)

%% Create Slices
plane_vec=[0,0,1]';
slice_sets={};
parfor count_slice=1:length(z_slices)
% for count_slice=1:length(z_slices)
        plane_point=[0,0,z_slices(count_slice)];
        current_vertex=geom_rot.vertices-plane_point;
        vert_dir=sign(current_vertex*plane_vec);
        elems_list=[];
        % the following lines determine all element with either a a vertext
        % on the plane, or part of the face cuts the plane
        for count_elems=1:size(geom_rot.faces,1)
                nodel=geom_rot.faces(count_elems,:);
                if vert_dir(nodel(1)) ~= vert_dir(nodel(2)) || ...
                                vert_dir(nodel(2)) ~= vert_dir(nodel(3))
                        elems_list=[elems_list;geom_rot.faces(count_elems,:)];
                elseif vert_dir(nodel(1)) ==0 &&...
                                vert_dir(nodel(2)) ==0 &&...
                                vert_dir(nodel(3)) ==0
                        elems_list=[elems_list;geom_rot.faces(count_elems,:)];
                end
        end
        % have a list of all elements that cut the current plane
        points_vec=[];
        segments_vec={};
        counter=1;
        for count_elems_plane=1:size(elems_list,1)
                nodel=elems_list(count_elems_plane,:);
                temp_vertices=geom_rot.vertices(nodel,:);
                current_vert_dir=vert_dir(nodel);
                 % the following statements check if the vertex is on the
                 % plane
%                 if current_vert_dir(1)==0
%                         points_vec=[points_vec;temp_vertices(1,:)];
%                 end
%                 if current_vert_dir(2)==0
%                         points_vec=[points_vec;temp_vertices(2,:)];
%                 end
%                 if current_vert_dir(3)==0
%                         points_vec=[points_vec;temp_vertices(3,:)];
%                 end
                % the following lines determine where the edge interesects
                % the plane
                current_pt=[];
                if current_vert_dir(1) ~= current_vert_dir(2)
                        intersect_point_x=interp1([temp_vertices(1,3),temp_vertices(2,3)],...
                                [temp_vertices(1,1),temp_vertices(2,1)],z_slices(count_slice));
                        intersect_point_y=interp1([temp_vertices(1,3),temp_vertices(2,3)],...
                                [temp_vertices(1,2),temp_vertices(2,2)],z_slices(count_slice));
                        intersect_point=[intersect_point_x,intersect_point_y,z_slices(count_slice)];
                        points_vec=[points_vec;intersect_point];
                        current_pt=[current_pt;intersect_point];
                end
                if current_vert_dir(1) ~= current_vert_dir(3)
                        intersect_point_x=interp1([temp_vertices(1,3),temp_vertices(3,3)],...
                                [temp_vertices(1,1),temp_vertices(3,1)],z_slices(count_slice));
                        intersect_point_y=interp1([temp_vertices(1,3),temp_vertices(3,3)],...
                                [temp_vertices(1,2),temp_vertices(3,2)],z_slices(count_slice));
                        intersect_point=[intersect_point_x,intersect_point_y,z_slices(count_slice)];
                        points_vec=[points_vec;intersect_point];
                        current_pt=[current_pt;intersect_point];
                end
                if current_vert_dir(2) ~= current_vert_dir(3)
                        intersect_point_x=interp1([temp_vertices(2,3),temp_vertices(3,3)],...
                                [temp_vertices(2,1),temp_vertices(3,1)],z_slices(count_slice));
                        intersect_point_y=interp1([temp_vertices(2,3),temp_vertices(3,3)],...
                                [temp_vertices(2,2),temp_vertices(3,2)],z_slices(count_slice));
                        intersect_point=[intersect_point_x,intersect_point_y,z_slices(count_slice)];
                        points_vec=[points_vec;intersect_point];
                        current_pt=[current_pt;intersect_point];
                end
                if ~isempty(current_pt)
                        segments_vec{counter}=current_pt;
                        counter=counter+1;
                end
        end
%         if ~isempty(points_vec)
%                 clf
%                 plot(points_vec(:,1),points_vec(:,2),'ro');
%                 xlim([min_bound_vertices(1), max_bound_vertices(1)]);
%                 ylim([min_bound_vertices(2), max_bound_vertices(2)]);
%                 axis equal
%                 hold on
%                 current_slice=segments_vec;
%                 for count_segment=1:length(current_slice)
%                         line_segment=current_slice{count_segment};
%                         plot(line_segment(:,1),line_segment(:,2),'k');
%                         hold on
%                 end
%                 pause(.001);
%         end
        slice_sets{count_slice}=segments_vec;
end
%% create boundaries plots
% for count_slice=185%1:length(slice_sets)
%         
%         clf
%         current_slice=slice_sets{count_slice};
%         for count_segment=1:length(current_slice)
%                 line_segment=current_slice{count_segment};
%                 plot(line_segment(:,1),line_segment(:,2),'k');
% 
%                 hold on
%         end
%         xlim([min_bound_vertices(1),max_bound_vertices(1)]);
%         ylim([min_bound_vertices(2),max_bound_vertices(2)]);
%         axis equal
%         pause(.001)
% end

%% determine segments for slice
for count_current_slice=1:length(slice_sets)
        current_slice=slice_sets{count_current_slice};
        
        plane_vertex=[];
        plane_segment_ids=[];
        slice_empty=0;
        if ~isempty(current_slice)
                vertex1=current_slice{1}(1,:);
                vertex2=current_slice{1}(2,:);
                plane_vertex=[plane_vertex;vertex1;vertex2];
                plane_segment_ids=[1,2];
                for count_segment=2:length(current_slice)
                        vertex1=current_slice{count_segment}(1,:);
                        vertex2=current_slice{count_segment}(2,:);
                        plane_index=[];
                        if ~ismember(vertex1,plane_vertex,'rows')
                                plane_vertex=[plane_vertex;vertex1];
                                plane_index=[plane_index,size(plane_vertex,1)];
                        else
                                [~,temp_index]=ismember(vertex1,plane_vertex,'rows');
                                plane_index=[plane_index,temp_index];
                        end
                        if ~ismember(vertex2,plane_vertex,'rows')
                                plane_vertex=[plane_vertex;vertex2];
                                plane_index=[plane_index,size(plane_vertex,1)];
                        else
                                [~,temp_index]=ismember(vertex2,plane_vertex,'rows');
                                plane_index=[plane_index,temp_index];
                        end
                        plane_segment_ids(count_segment,:)=plane_index;
                end
        else
                slice_empty=1;
        end
        %% Create loop segments
        if ~slice_empty
                [loop_segments,~,~]=getClosedLoops(plane_segment_ids);
                if exist('vertex1','var')
                        if logical(min(vertex1==vertex2))
                                slice_empty=1;
                        end
                end
        end
        %% Create Loop Node List
        if ~slice_empty
                loop_nodes_segments={};
                for counti=1:length(loop_segments)
                        loop_nodes_segments{counti}=plane_vertex(loop_segments{counti}',:);
                end
                slices_loops{count_current_slice}=loop_nodes_segments;
                %% Verify segment loops
%                 hold on
%                 for count_loop=1:length(loop_nodes_segments)
%                         plot(loop_nodes_segments{count_loop}(:,1),loop_nodes_segments{count_loop}(:,2),'o');
%                         axis equal
%                         hold on
%                 end

                %% create polygons
                clear slice_polygons
                for counti_loops=1:length(loop_nodes_segments)
                      slice_polygons(counti_loops)=polyshape(loop_nodes_segments{counti_loops}(:,1),loop_nodes_segments{counti_loops}(:,2));
                end
                slices_polygons{count_current_slice}=slice_polygons;
                %% Determine Polygons inside
                polygon_inside_mat=zeros(length(slice_polygons));
                for counti_loops=1:length(slice_polygons)
                        for countj_loops=1:length(slice_polygons)
                                ps1=slice_polygons(counti_loops);
                                if counti_loops ~= countj_loops
                                        ps2=slice_polygons(countj_loops);
                                        isInside = abs(area(ps2)-area(subtract(ps2,ps1))-area(ps1))<(area(ps1)*1e-06);
                                        if isInside
                                                polygon_inside_mat(counti_loops,countj_loops)=isInside;
                                        end
                                end
                        end
                end
                polygon_inside_mat2=polygon_inside_mat;
                max_iters=100;
                counter=1;

                tree_array=zeros(length(polygon_inside_mat2),1);
                while sum(sum(polygon_inside_mat2))>0 && counter<max_iters
                        for count_poly=1:length(polygon_inside_mat2)
                                if sum(polygon_inside_mat2(count_poly,:))==1
                                        in_poly_id=find(polygon_inside_mat2(count_poly,:)==1);
                                        tree_array(count_poly)=in_poly_id;
                                        polygon_inside_mat2(count_poly,:)=0;
                                        for count_polyj=1:length(polygon_inside_mat2)
                                                if polygon_inside_mat2(count_polyj,count_poly) ==1
                                                        polygon_inside_mat2(count_polyj,in_poly_id)=0;
                                                end
                                        end
                                end
                        end
                        counter=counter+1;
                end
%                 treefig=figure();
%                 treeplot(tree_array');
%                 [tx,ty]=treelayout(tree_array);
%                 text(tx+0.02,ty,num2cell([1:length(tree_array)]))

                %% Determine Polygon Type
                rank_array=tree_array;
                for count_poly=1:length(tree_array)
                        current_lev=tree_array(count_poly);
                        rank_lev=0;
                        while current_lev ~= 0
                                current_lev=tree_array(current_lev);
                                rank_lev=rank_lev+1;
                        end
                        rank_array(count_poly)=rank_lev;
                end

%                 poly_shading_array=mod(rank_array+1,2);

                %% Create Image
                if ~isempty(slice_polygons)
                        treemat=zeros(length(tree_array));
                        for counti=1:length(tree_array)
                                if tree_array(counti)~=0
                                        treemat(counti,tree_array(counti))=1;
                                end
                        end
                        slice_bin_vol = mpoly2mask({slice_polygons.Vertices}, x_slices, (y_slices),treemat);
                        slices_bin_vol(:,:,count_current_slice)=slice_bin_vol;
                end
        else
                slices_loops{count_current_slice}=[];
                slices_polygons{count_current_slice}=[];
        end
end
%% plot slices
% disp('total Image binarization time');
% 
% slice_figure=figure('WindowState','Maximized');
% for counti=1:size(slices_bin_vol,3)
%       
%       slice_plot=imshow(slices_bin_vol(:,:,counti),'Parent',gca);
% %       slice_figure.set('WindowState','Maximized');
%       pause(.001);
% end

%% Save current set
% save([filename,'sliced.mat'],'slices_bin_vol','slices_polygons','slices_loops',...
%         'pixel_bounds','edge_bounds','x_spacing','y_spacing','z_spacing',...
%         'slice_direction_z','geom','x_slices','y_slices','z_slices','slice_dir_RotMat',...
%         'vmin_orig','vmax_orig','vmin_rot','vmax_rot')