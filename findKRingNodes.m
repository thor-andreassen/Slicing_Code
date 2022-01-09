function [k_node_list,ring_nodes_layered]=findKRingNodes(faces,center_node,k_rings)
        k=1;
        k_node_list=[center_node];
        ring_nodes_layered{1}=[center_node];
        while k<=k_rings
                current_node_layer=[];
                for node_num=1:length(ring_nodes_layered{k})
                        current_node=ring_nodes_layered{k}(node_num);
                        [faces_with_node,~]=find(faces==current_node);
                        current_nodes=zeros(length(faces_with_node),3);
                        for count_face=1:length(faces_with_node)
                                current_nodes(count_face,:)=faces(faces_with_node(count_face),:);
                        end
                        current_nodes_long=reshape(current_nodes,[length(faces_with_node)*3,1]);
                        current_nodes_long=unique(current_nodes_long);
                        current_id=find(current_nodes_long==node_num);
                        current_nodes_long(current_id)=[];
                        current_node_layer=[current_node_layer;current_nodes_long];
                end
                for count_node=1:length(k_node_list)
                        current_id=find(current_node_layer==k_node_list(count_node));
                        if ~isempty(current_id)
                                current_node_layer(current_id)=[];
                        end
                end
                k_node_list=[k_node_list;current_node_layer];
                ring_nodes_layered{k+1}=current_node_layer;
                k=k+1;
        end
end