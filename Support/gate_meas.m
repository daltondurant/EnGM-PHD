% Gate Measurements using the Squared Mahalanobis Distance
function z_gate = gate_meas(cfig,model,z,m,P)
%
    valid_idx = [];
    zlength = size(z,2); 
    if zlength == 0 
        z_gate= []; 
        return; 
    end
    plength = size(m,2);
    
    for jj=1:plength
        Hj = cfig.H(model,m(:,jj));
        Sj= Hj*P(:,:,jj)*Hj' + model.R; 
        Sj= (Sj+ Sj')/2; 
        Vs= chol(Sj); 
        inv_sqrt_Sj= inv(Vs);
        nu= z - repmat(cfig.h(model,m(:,jj),zeros(size(model.D,2),1)),[1 zlength]);
        dist = sum((inv_sqrt_Sj'*nu).^2); % squared Mahalanobis Distance 
        valid_idx = union(valid_idx,find( dist < model.gamma_G ));
    end
    z_gate = z(:,valid_idx);
%
end