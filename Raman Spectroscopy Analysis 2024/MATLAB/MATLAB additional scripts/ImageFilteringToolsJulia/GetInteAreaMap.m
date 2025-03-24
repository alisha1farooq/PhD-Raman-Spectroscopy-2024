function [ChemicalMap2D,ChemicalMap3D] = GetInteAreaMap(wdf)

clc
count = wdf.Count();
    ReshapeSize  = wdf.GetWMapblock().numPoints;
    ReshapeSize(3) = [];
    ReshapeSize  = ReshapeSize';
Spectra = zeros(count,1);
l=1;
wdf.StartChunkwiseReading();
while (wdf.AreMoreChunks());
    z = wdf.GetNextDataChunk();
    Size = size(z);
    numberReadIn = Size(1);
    Spectra(l:l+numberReadIn-1) = trapz(z,2);%sum(z2);%cumtrapz(z,2);
    l = l + numberReadIn;
end

ChemicalMap2D=Spectra;
ChemicalMap3D=reshape(ChemicalMap2D,ReshapeSize(1),ReshapeSize(2),1);

end