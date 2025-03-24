function RaMP = WriteDatatoRaMP(XList,Map,ReshapeSize,mask)
%
RaMP=[];
RaMP.xAxisData=XList;
RaMP.dataDims=[ReshapeSize(1),ReshapeSize(2),size(XList,2)];
RaMP.matrix=Map;
RaMP.toRemove=zeros(ReshapeSize(1),ReshapeSize(2),1);
RaMP.toRemove(mask')=1;
RaMP.IN=[];
RaMP.xj=XList;
end