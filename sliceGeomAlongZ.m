function [slice_sets,slices_loops,slices_polygons,slices_bin_vol]=...
sliceGeomAlongZ(geom,x_slices, y_slices,z_slices,slice_bin_vol,slices_bin_vol)
        %% Main Slicing Function to slice a geometry according to the slice parameters and create a binary mesh
        % Created by Thor Andreassen
        % 1/10/22
        %
        %
        % This function uses a series of steps to slice a geometry into a series of binary matrices
        % The function takes in a geometry that is assumed to be aligned
        % with the desired slice direction as the "z-direction" of the
        % mesh, and then takes in parameters defining the mesh slices
        % locations. The result is a set of slice geometries and contour
        % polygons as well as a matrix of binary values for regions
        % "inside" and "outside" the slice.
        %
        %
        % Inputs:
                % geom: a structure defining the current mesh surface to be
                        % sliced as the following:
                                % geom.faces = (n x 3) connectivity list of
                                        % triangulated mesh curfaces defining the
                                        % corner nodes of each face as a row
                                % geom.vertices = (m x 3) matrix of x, y, z
                                        % coordinates of all vertices/nodes used to
                                        % define the mesh surface. Each row
                                        % corresponds to the node number of the
                                        % mesh defined in "faces"
                % x_slices: a vector containing all of the x coordinates to
                        % turn into the x- spacing pixels of the final image.
                % y_slices: a vector containing all of the y coordinates to
                        % turn into the y- spacing pixels of the final image.
                % z_slices: a vector containing all of the z coordinates
                        % that define the resulting slices of the algorithm
                % slice_bin_vol: a previously define binary logic matrix
                        % representing a single slice of the resulting slices
                        % the size of this matrix is (length(x_slices) x length(y_slices))
                        % (optional)
                % slices_bin_vol a previously define binary logic matrix
                        % representing all of the slice values of the resulting slices
                        % the size of this matrix is (length(x_slices) x length(y_slices) x length(z_slices))
                        % (optional)
        %
        %
        %
        % Outputs:
                % slice_sets: A cell array containing all of the points
                        % (x,y,z) that intersect with the current z slice. 
                        % The size of this matrix is the same as the length of
                        % z_slices
                % slices_loops: A cell array containing the loop sets for
                        % each slice. Where each cell contains a cell array with
                        % all of the closed loops separated. % The size of this matrix
                        % is the same as the length of z_slices
                % slices_polygons: A cell array containing all of the
                        % polygons for each slice. % The size of this matrix is the
                        % same as the length of z_slices
                % slices_bin_vol: The resulting binary matrix of all of the
                        % slices of the geometry. The size is (length(x_slices)+1 x
                        % length (y_slices)+1 x length (z_slices)). the
                        % matrix contains a 1 wherever the voxel is
                        % contained in a solid part of the original
                        % geometry, and a 0 wherever there is no geometry
                        % present.
                        
                        % the following line is used to initialize the
                        % outputs if they are not previously passed in from
                        % an additive matrix from the calling function.
        if nargin < 5
                x_spacing=x_slices(2)-x_slices(1);
                y_spacing=y_slices(2)-y_slices(1);
                z_spacing=z_slices(2)-z_slices(1);
                max_bound_vertices=[x_slices(end),y_slices(end),z_slices(end)];
                min_bound_vertices=[x_slices(1),y_slices(1),z_slices(1)];
                
                volume_pix_num=ceil((max_bound_vertices-min_bound_vertices)./[x_spacing,y_spacing,z_spacing]);
                total_pix_num=volume_pix_num;
                slice_bin_vol=false(total_pix_num(2)+1,total_pix_num(1)+1);
                slices_bin_vol=logical(zeros(total_pix_num(2)+1,total_pix_num(1)+1,length(z_slices)));
        end
        %% Create Slices
        % the following section is the main portion of the code to slice
        % the geometry along z_Sliced planes by determining the
        % intersection between the faces of the mesh and a plane made at
        % the current slice.
        use_parallel=1;
        [slice_sets]=getSliceSetsInZ(geom,z_slices,use_parallel);
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
        % the following section is the main piece of the code that takes
        % the sliced coordinates across the geometry, and converts it to
        % closed regions that are then binarized accordingly.
        for count_current_slice=1:length(slice_sets)
                
                current_slice=slice_sets{count_current_slice};
                
                % the following section of code is used to take the set of
                % nodes from the line edges of the slice, and convert them
                % to a series of "new nodes" with edges defined as
                % connecting the set of node numbers together, instead of
                % as 2 sets of (x, y, z) coordinates. The loop goes through
                % each node vertex in the list of line segments, and if the
                % node coordinates already exist in the assigned list of
                % nodes, it assigns the node number to it, and continues.
                % If it does not already exist, it creates a new node as
                % the next integer in the list, and adds the connection.
                % This is done for both ends of each line segment, and then
                % repeats for each line segment.
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
                % the following section of code works on the list of node
                % connections, and not the vertex coordinates themselves,
                % and determines the closed set of loops defined by those
                % nodes. The loop goes through each slice and calls the
                % separate function getClosedLoops to find the closed
                % sections.
                if ~slice_empty
                        plane_segment_ids=unique(plane_segment_ids,'rows');
                        for count_ids=size(plane_segment_ids,1):-1:1
                                if plane_segment_ids(count_ids,1)==plane_segment_ids(count_ids,2)
                                        plane_segment_ids(count_ids,:)=[];
                                end
                        end
                        [loop_segments,closed_loops,~]=getClosedLoops(plane_segment_ids);
                        for count_loop=length(closed_loops):-1:1
                                if closed_loops(count_loop)==0
                                        loop_segments(count_loop)=[];
                                end
                        end
                        if exist('vertex1','var')
                                if logical(min(vertex1==vertex2))
                                        slice_empty=1;
                                end
                        end
                end
                %% Create Loop Node List
                % the following set of code takes the previously determined
                % set of closed loops from the nodes, and converts them to
                % closed sets of x,y,z vertices based on teh nodes and
                % their coordinates.
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
                        % the following lines create matlab polygon shapes
                        % based on the closed loop vertices, so as to allow
                        % easy checking for nodes inside and outside the
                        % polygons.
                        clear slice_polygons
                        for counti_loops=1:length(loop_nodes_segments)
                              slice_polygons(counti_loops)=polyshape(loop_nodes_segments{counti_loops}(:,1),loop_nodes_segments{counti_loops}(:,2));
                        end
                        if exist('slice_polygons','var')
                                for count_poly=length(slice_polygons):-1:1
                                        if isempty(slice_polygons(count_poly).Vertices)
                                                slice_polygons(count_poly)=[];
                                        end
                                end
                                slices_polygons{count_current_slice}=slice_polygons;
                        end
                        %% Determine Polygons inside
                        % the following loops determine for the current
                        % slice which of the closed loop polygons appears
                        % in which of the others. The loop goes through
                        % each of the polygons, (assuming they are
                        % different) and creates a matrix that describes if
                        % the polygon is completely contained within the
                        % other polygon. If it is, then there is a 1 in
                        % that space and if not, there is a 0. This matrix
                        % can then be used to determine whether a given
                        % polygon should be "filled" or whether it is an
                        % internal hole.
                        polygon_inside_mat=[];
                        if exist('slice_polygons','var')
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
                        end
                        polygon_inside_mat2=polygon_inside_mat;
                        max_iters=100;
                        counter=1;
                        
                        % the following takes the previously created
                        % matrix of "inside" polygons for the current
                        % slice, and breaks it apart into which polygons
                        % are direclty inside of which. In particular it
                        % then determine based on the number of polygons it
                        % is inside, the rank of the polygon and
                        % subsequently whether the polygon is "filled" or
                        % "unfilled".
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
                        % the following determine the rank/level of the
                        % given polygon which is used to determine the
                        % color based on whether it is odd or even.
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
                        % the following line calls a function that takes
                        % the set of polygons and the chosen slice
                        % dimensions to create a single binary image for
                        % the slice.
                        if exist('slice_polygons','var')
                                if ~isempty(slice_polygons)
                                        treemat=zeros(length(tree_array));
                                        for counti=1:length(tree_array)
                                                if tree_array(counti)~=0
                                                        treemat(counti,tree_array(counti))=1;
                                                end
                                        end
                                        try
                                                slice_bin_vol = mpoly2mask({slice_polygons.Vertices}, x_slices, (y_slices),treemat);
                                        catch
                                                 disp('polymaskfailed');
                                        end
                                        slices_bin_vol(:,:,count_current_slice)=slice_bin_vol;
                                end
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
end