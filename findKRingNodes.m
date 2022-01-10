function [k_node_list,ring_nodes_layered]=findKRingNodes(faces,center_node,k_rings)
        %% Function to list N connected Nodes to a given node
        % Created by Thor Andreassen
        % 1/10/22
        % The function counts through a list of faces, and determines every node
        % that is connected to the chosen center node. Then the nodes creates a
        % ring of the nodes that are connected to the previous ring of nodes based
        % on the nodes attache to the face with the previous ring of nodes. The
        % ring number is equivalent to determining the number of edges between the
        % enter node, and each of the nodes in subsequent rings.
        %
        %
        % Inputs:
                % faces: a list of (n x 3) connected triangluar meshes as the node
                        % numbers of the nodes on the corner of each face
                % center_node: The node number corresponding to the center node of
                        % the chosen ring set of nearest ring neighbors.
                % k_rings: The number of  rings to find.
                        % K=1 contains only the original center node.
                        % K=2 is the first layer of connected nodes, based
                                % on faces attached to the center node.
                        % K=... additional layers
        %
        %
        % Outputs:
                % k_node_list: The set of all nodes that are within K rings of the
                        % center node.
                % ring_nodes_layered: A cell array containing all of the nodes in
                        % each ring separated into different rings. 
                k=1;
                k_node_list=[center_node];
                ring_nodes_layered{1}=[center_node];
                % loop goes through each of the rings
                while k<=k_rings
                        current_node_layer=[];
                        % the loop goes through all of the nodes that defined the
                        % previous ring.
                        for node_num=1:length(ring_nodes_layered{k})
                                % the loop chooses the current node form the list
                                % of nodes in the previous ring
                                current_node=ring_nodes_layered{k}(node_num);

                                % the following finds all of the faces that contain
                                % the current node.
                                [faces_with_node,~]=find(faces==current_node);

                                % the following loop goes through all of the faces
                                % that include that node, and add the nodes to a
                                % list.
                                current_nodes=zeros(length(faces_with_node),3);
                                for count_face=1:length(faces_with_node)
                                        current_nodes(count_face,:)=faces(faces_with_node(count_face),:);
                                end

                                % the following lines remove any duplicate nodes
                                % from the current ring, if they occur on multiple
                                % faces.
                                current_nodes_long=reshape(current_nodes,[length(faces_with_node)*3,1]);
                                current_nodes_long=unique(current_nodes_long);
                                current_id=find(current_nodes_long==node_num);
                                current_nodes_long(current_id)=[];
                                current_node_layer=[current_node_layer;current_nodes_long];
                        end
                        % the following lines remove any nodes that exist in
                        % previous rings.
                        for count_node=1:length(k_node_list)
                                current_id=find(current_node_layer==k_node_list(count_node));
                                if ~isempty(current_id)
                                        current_node_layer(current_id)=[];
                                end
                        end
                        % the following lines sets the full node list to the current
                        % ring and add thats that to the list of rings, as well as
                        % the total set of nodes within the number of rings.
                        k_node_list=[k_node_list;current_node_layer];
                        ring_nodes_layered{k+1}=current_node_layer;
                        k=k+1;
                end
        end