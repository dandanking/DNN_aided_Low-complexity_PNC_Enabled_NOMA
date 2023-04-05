
function phyBridgeMap = phy_gen_bridge_map(phyRawMap)

phyBridgeMap = zeros(length(phyRawMap),3);
A_INDEX=1; B_INDEX=2; X_INDEX=3;
for jj=1:length(phyRawMap)
    okA = phyRawMap(jj,A_INDEX);
    okB = phyRawMap(jj,B_INDEX);
    okX = phyRawMap(jj,X_INDEX);
    if (okA&&~okB&&okX) || (~okA&&okB&&okX)
        okA=1; okB=1;
    end
    phyBridgeMap(jj,[A_INDEX,B_INDEX,X_INDEX]) = [okA,okB,okX];
end

end