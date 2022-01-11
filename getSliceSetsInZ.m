function [slice_sets]=getSliceSetsInZ(geom,z_slices,use_parallel)
%% Function to slice the geometry in the z direction
        % created by Thor Andreassen
        % 1/10/22
        %
        % This function takes in a triangulated geometry and a set of
        % z-slices as z coordinates and determine the set of slice sets
        % that include all the points that intersect each of the z-slice
        % planes.
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
                % z_slices: a vector containing all of the z coordinates
                        % that define the resulting slices of the algorithm
                % use_parallel: the value is 1, if a parallel loop is to be
                % used, and 0 otherwise.
        % Outputs:
                %slice_sets: a cell array that has the same number of rows
                %as the z_slices and contains all the intersection line
                %segments at the current slice.
        plane_vec=[0,0,1]';
        slice_sets={};
        if nargin <=2
                use_parallel=0;
        end
        min_z=min(geom.vertices);
        min_z=min_z(end);
        
        max_z=max(geom.vertices);
        max_z=max_z(end);
        if use_parallel==1
                parfor count_slice=1:length(z_slices)
                        if z_slices(count_slice)>=min_z && z_slices(count_slice)<=max_z
                                % for count_slice=1:length(z_slices)

                                % the following lines define the point on the corner of the
                                % plane of the current intersection slice, as x=0, y=0 and
                                % z = the current slice coordinate.
                                plane_point=[0,0,z_slices(count_slice)];

                                % the following line determine the cartesian distance to
                                % the point in the plane by subtracting all the points of
                                % the geometry by the plane point.
                                current_vertex=geom.vertices-plane_point;

                                % the following line determines the sign of the distance to
                                % each of the point in the vertex. to the z direction.
                                vert_dir=sign(current_vertex*plane_vec);
                                elems_list=[];
                                % the following lines determine all element with either a a vertext
                                % on the plane, or part of the face cuts the plane
                                for count_elems=1:size(geom.faces,1)

                                        % the following lines go through each of the faces,
                                        % and check the signs of the distances at each of
                                        % the nodes. If at least one of the nodes on the
                                        % face has a positive distance and the other is
                                        % negative, then the face intersects the plane as
                                        % one point is above and one point is below. If so,
                                        % it adds the face to a list of faces that
                                        % intersect the plane.
                                        nodel=geom.faces(count_elems,:);
                                        if vert_dir(nodel(1)) ~= vert_dir(nodel(2)) || ...
                                                        vert_dir(nodel(2)) ~= vert_dir(nodel(3))
                                                elems_list=[elems_list;geom.faces(count_elems,:)];

                                                % the following portion determines if one of the
                                                % nodes is on the plane by having a zero distance.
                                                % if so it adds the face to the list of the faces
                                                % interseting the plane.
                                        elseif vert_dir(nodel(1)) ==0 &&...
                                                        vert_dir(nodel(2)) ==0 &&...
                                                        vert_dir(nodel(3)) ==0
                                                elems_list=[elems_list;geom.faces(count_elems,:)];
                                        end
                                end
                                % have a list of all elements that cut the current plane
                                segments_vec={};
                                counter=1;
                                % the following loop goes through the list of all the
                                % elements that intersect the plane and determine the
                                % intersection points.
                                for count_elems_plane=1:size(elems_list,1)
                                        % the following line determine the nodes of the
                                        % current intersection face.
                                        nodel=elems_list(count_elems_plane,:);
                                        temp_vertices=geom.vertices(nodel,:);
                                        current_vert_dir=vert_dir(nodel);
                                        % the following lines determine where the edge interesects
                                        % the plane
                                        current_pt=[];
                                        % the following lines determine which of the pairs
                                        % of nodes of the triangulated mesh cross the
                                        % plane. For each pair of points where one point is
                                        % above the plane, and another point is below the
                                        % plane, the point that intersects the plane is
                                        % determined by the interpolationg of the points to
                                        % the z point on the plane. The point is then
                                        % calculated as the x,y,z location of the
                                        % intersection, and then added to a list of points
                                        % for the current slice, as well as a total list
                                        % across all slices.
                                        if current_vert_dir(1) ~= current_vert_dir(2)
                                                intersect_point_x=interp1([temp_vertices(1,3),temp_vertices(2,3)],...
                                                        [temp_vertices(1,1),temp_vertices(2,1)],z_slices(count_slice));
                                                intersect_point_y=interp1([temp_vertices(1,3),temp_vertices(2,3)],...
                                                        [temp_vertices(1,2),temp_vertices(2,2)],z_slices(count_slice));
                                                intersect_point=[intersect_point_x,intersect_point_y,z_slices(count_slice)];
                                                current_pt=[current_pt;intersect_point];
                                        end
                                        if current_vert_dir(1) ~= current_vert_dir(3)
                                                intersect_point_x=interp1([temp_vertices(1,3),temp_vertices(3,3)],...
                                                        [temp_vertices(1,1),temp_vertices(3,1)],z_slices(count_slice));
                                                intersect_point_y=interp1([temp_vertices(1,3),temp_vertices(3,3)],...
                                                        [temp_vertices(1,2),temp_vertices(3,2)],z_slices(count_slice));
                                                intersect_point=[intersect_point_x,intersect_point_y,z_slices(count_slice)];
                                                current_pt=[current_pt;intersect_point];
                                        end
                                        if current_vert_dir(2) ~= current_vert_dir(3)
                                                intersect_point_x=interp1([temp_vertices(2,3),temp_vertices(3,3)],...
                                                        [temp_vertices(2,1),temp_vertices(3,1)],z_slices(count_slice));
                                                intersect_point_y=interp1([temp_vertices(2,3),temp_vertices(3,3)],...
                                                        [temp_vertices(2,2),temp_vertices(3,2)],z_slices(count_slice));
                                                intersect_point=[intersect_point_x,intersect_point_y,z_slices(count_slice)];
                                                current_pt=[current_pt;intersect_point];
                                        end
                                        % the following line checks to see if the plane
                                        % intersected and adds the points to a list of the
                                        % line segments present.
                                        if ~isempty(current_pt)
                                                segments_vec{counter}=current_pt;
                                                counter=counter+1;
                                        end
                                end
                                slice_sets{count_slice}=segments_vec;
                        else
                                slice_sets{count_slice}=[];
                        end
                end
        else
                for count_slice=1:length(z_slices)
                        if z_slices(count_slice)>=min_z && z_slices(count_slice)<=max_z
                                % for count_slice=1:length(z_slices)

                                % the following lines define the point on the corner of the
                                % plane of the current intersection slice, as x=0, y=0 and
                                % z = the current slice coordinate.
                                plane_point=[0,0,z_slices(count_slice)];

                                % the following line determine the cartesian distance to
                                % the point in the plane by subtracting all the points of
                                % the geometry by the plane point.
                                current_vertex=geom.vertices-plane_point;

                                % the following line determines the sign of the distance to
                                % each of the point in the vertex. to the z direction.
                                vert_dir=sign(current_vertex*plane_vec);
                                elems_list=[];
                                % the following lines determine all element with either a a vertext
                                % on the plane, or part of the face cuts the plane
                                for count_elems=1:size(geom.faces,1)

                                        % the following lines go through each of the faces,
                                        % and check the signs of the distances at each of
                                        % the nodes. If at least one of the nodes on the
                                        % face has a positive distance and the other is
                                        % negative, then the face intersects the plane as
                                        % one point is above and one point is below. If so,
                                        % it adds the face to a list of faces that
                                        % intersect the plane.
                                        nodel=geom.faces(count_elems,:);
                                        if vert_dir(nodel(1)) ~= vert_dir(nodel(2)) || ...
                                                        vert_dir(nodel(2)) ~= vert_dir(nodel(3))
                                                elems_list=[elems_list;geom.faces(count_elems,:)];

                                                % the following portion determines if one of the
                                                % nodes is on the plane by having a zero distance.
                                                % if so it adds the face to the list of the faces
                                                % interseting the plane.
                                        elseif vert_dir(nodel(1)) ==0 &&...
                                                        vert_dir(nodel(2)) ==0 &&...
                                                        vert_dir(nodel(3)) ==0
                                                elems_list=[elems_list;geom.faces(count_elems,:)];
                                        end
                                end
                                % have a list of all elements that cut the current plane
                                segments_vec={};
                                counter=1;
                                % the following loop goes through the list of all the
                                % elements that intersect the plane and determine the
                                % intersection points.
                                for count_elems_plane=1:size(elems_list,1)
                                        % the following line determine the nodes of the
                                        % current intersection face.
                                        nodel=elems_list(count_elems_plane,:);
                                        temp_vertices=geom.vertices(nodel,:);
                                        current_vert_dir=vert_dir(nodel);
                                        % the following lines determine where the edge interesects
                                        % the plane
                                        current_pt=[];
                                        % the following lines determine which of the pairs
                                        % of nodes of the triangulated mesh cross the
                                        % plane. For each pair of points where one point is
                                        % above the plane, and another point is below the
                                        % plane, the point that intersects the plane is
                                        % determined by the interpolationg of the points to
                                        % the z point on the plane. The point is then
                                        % calculated as the x,y,z location of the
                                        % intersection, and then added to a list of points
                                        % for the current slice, as well as a total list
                                        % across all slices.
                                        if current_vert_dir(1) ~= current_vert_dir(2)
                                                intersect_point_x=interp1([temp_vertices(1,3),temp_vertices(2,3)],...
                                                        [temp_vertices(1,1),temp_vertices(2,1)],z_slices(count_slice));
                                                intersect_point_y=interp1([temp_vertices(1,3),temp_vertices(2,3)],...
                                                        [temp_vertices(1,2),temp_vertices(2,2)],z_slices(count_slice));
                                                intersect_point=[intersect_point_x,intersect_point_y,z_slices(count_slice)];
                                                current_pt=[current_pt;intersect_point];
                                        end
                                        if current_vert_dir(1) ~= current_vert_dir(3)
                                                intersect_point_x=interp1([temp_vertices(1,3),temp_vertices(3,3)],...
                                                        [temp_vertices(1,1),temp_vertices(3,1)],z_slices(count_slice));
                                                intersect_point_y=interp1([temp_vertices(1,3),temp_vertices(3,3)],...
                                                        [temp_vertices(1,2),temp_vertices(3,2)],z_slices(count_slice));
                                                intersect_point=[intersect_point_x,intersect_point_y,z_slices(count_slice)];
                                                current_pt=[current_pt;intersect_point];
                                        end
                                        if current_vert_dir(2) ~= current_vert_dir(3)
                                                intersect_point_x=interp1([temp_vertices(2,3),temp_vertices(3,3)],...
                                                        [temp_vertices(2,1),temp_vertices(3,1)],z_slices(count_slice));
                                                intersect_point_y=interp1([temp_vertices(2,3),temp_vertices(3,3)],...
                                                        [temp_vertices(2,2),temp_vertices(3,2)],z_slices(count_slice));
                                                intersect_point=[intersect_point_x,intersect_point_y,z_slices(count_slice)];
                                                current_pt=[current_pt;intersect_point];
                                        end
                                        % the following line checks to see if the plane
                                        % intersected and adds the points to a list of the
                                        % line segments present.
                                        if ~isempty(current_pt)
                                                segments_vec{counter}=current_pt;
                                                counter=counter+1;
                                        end
                                end
                                slice_sets{count_slice}=segments_vec;
                        else
                                slice_sets{count_slice}=[];
                        end
                end
        end
end