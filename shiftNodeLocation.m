function [new_node_loc,disp]=shiftNodeLocation(faces,nodes,num_node,...
        near_nodes_id,ring_nodes_layered)
        % Created by Thor Andreassen
        % 1/10/22
        % This code implements an algorithm by Wang et al 2011 called:
        % "Quality Mesh Smoothing via Local Surface Fitting and Optimum"
        %
        % The function uses a set of faces, nodes and connections between
        % nodes to approximate the surface at each node by fitting a
        % quadric surface to the nodes. Then the modely projects these
        % poitns onto a plane for the node vertex based on an approximate
        % by its attached faces, and then determines the "ideal" node
        % location in the projected surface by solving the maximum
        % inscribed triangle problem. The new projected location is then
        % calculated relative to the quadric surface fitted for the
        % original nodes to calculate the new optimal location for the
        % node. This is then repeated across all nodes to find the correct
        % node location.
        %
        % Inputs:
                % faces: (n x 3) matrix of connections of nodes defined by
                        % corner node indices for a triangulated mesh
                % nodes: (mx3) matrix of node coordinates (x, y, z) with each row
                        % corresponding to the assumed node number and
                        % corresponding to the nodes used in the "faces" input
                % node_num: a single number correpsonding to the node to
                        % calculate the normal of. this number should
                        % corresponding to the row of the node in the nodes
                        % input
                % near_nodes_id: An array that contains the node numbers
                        % of attached nodes to the current surface
                % ring_nodes_layered: A cell array that contains the
                        % adjacent rings attached to the current node.
        % Outputs:
                % new_node_loc: (x, y, z) coordinates of the new "smoothed"
                        % location of the chosen node.
                % disp: the amount that the node moved from its original
                        % location to is new smoothed location
        near_nodes=nodes(near_nodes_id,:);
%         [near_nodes_id,ring_nodes_layered]=findKRingNodes(faces,num_node,num_neighbors);

        vecz=findVertexNormalFromMesh(faces,nodes,num_node);
        temp_vecx=near_nodes(2,:)-near_nodes(1,:);
        temp_vecx=temp_vecx/norm(temp_vecx);
        vecy=cross(vecz,temp_vecx)/norm(cross(vecz,temp_vecx));
        vecx=cross(vecy,vecz)/norm(cross(vecy,vecz));

        rotmat=[vecx', vecy', vecz'];
        transvec=[nodes(num_node,:)'];

        LFCS_to_GCS_TransMat=eye(4);
        LFCS_to_GCS_TransMat(1:3,1:3)=rotmat;
        LFCS_to_GCS_TransMat(1:3,4)=transvec;
        GCS_to_LFCS_TransMat=inv(LFCS_to_GCS_TransMat);

        near_nodes_surf_project=GCS_to_LFCS_TransMat*[near_nodes';ones(1,size(near_nodes,1))];
        near_nodes_surf_project=near_nodes_surf_project(1:3,:)';


        nearest_ring_id=ring_nodes_layered{2};
        nearest_ring_nodes=nodes(nearest_ring_id,:);
        near_nodes_polygon_project=GCS_to_LFCS_TransMat*[nearest_ring_nodes';ones(1,size(nearest_ring_nodes,1))];
        near_nodes_polygon_project=near_nodes_polygon_project(1:3,:)';

        temp_res=max_inscribed_circle(near_nodes_polygon_project(:,1),near_nodes_polygon_project(:,2));
        new_node_loc_proj(1)=temp_res.xc;
        new_node_loc_proj(2)=temp_res.yc;
        [quadric_coefs,~]=fitQuadricSurf(near_nodes_surf_project(:,1),near_nodes_surf_project(:,2),near_nodes_surf_project(:,3));

        quadric_vec=[new_node_loc_proj(1).^2,...
                new_node_loc_proj(1).*new_node_loc_proj(2),...
                new_node_loc_proj(2).^2,...
                new_node_loc_proj(1),...
                new_node_loc_proj(2),...
                1];
        new_node_loc_proj(3)=quadric_vec*quadric_coefs;
        new_node_loc_GCS=LFCS_to_GCS_TransMat*[new_node_loc_proj';ones(1,size(new_node_loc_proj,1))];
        new_node_loc_GCS=new_node_loc_GCS(1:3)';
        new_node_loc=new_node_loc_GCS;
        disp=norm(new_node_loc-near_nodes(1,:));