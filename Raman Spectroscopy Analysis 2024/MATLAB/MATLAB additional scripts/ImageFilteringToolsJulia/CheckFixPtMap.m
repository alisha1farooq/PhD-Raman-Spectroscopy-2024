function [newSpectra,Map,ReshapeSize,mask]=CheckFixPtMap(Spec,fullFilename,mask)
%
tempwdf=WdfReader(fullFilename);
ReshapeSize  = tempwdf.GetWMapblock().numPoints;
ReshapeSize(3) = [];
ReshapeSize  = ReshapeSize';
XList = tempwdf.GetXList;
Count=tempwdf.Count;
X = tempwdf.GetOriginListValues(WiREDataType.SpatialX, 1, Count);
Y = tempwdf.GetOriginListValues(WiREDataType.SpatialY, 1, Count);
if tempwdf.MeasurementType==3&&ReshapeSize(1)==1&&ReshapeSize(2)==1
       a=floor(sqrt(Count));
       b=floor(Count/a);
       ReshapeSize(1)=a;
       ReshapeSize(2)=b;
       TotalX=a;
       TotalY=b;
       newSpectra=zeros(ReshapeSize(1)*ReshapeSize(2),size(Spec,2));
       newSpectra(1:Count,:)=Spec;
       addedmask=[Count+1:size(Spec,1)];
       mask=[mask;addedmask'];
else
    TotalX=unique(Y);TotalX=size(TotalX,1);
    TotalY=unique(X);TotalY=size(TotalY,1);
    newSpectra=Spec;
end

Map=reshape(newSpectra,ReshapeSize(1),ReshapeSize(2),size(Spec,2));     
tempwdf.Close();

end