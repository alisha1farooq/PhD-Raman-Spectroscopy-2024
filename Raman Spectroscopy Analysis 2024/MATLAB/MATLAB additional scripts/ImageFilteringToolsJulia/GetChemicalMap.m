function [Spectra,ChemicalMap2D,ChemicalMap3D] = GetChemicalMap(wdf,wavenumber)

clc

    XList=wdf.GetXList;
    count = wdf.Count();
    ReshapeSize  = wdf.GetWMapblock().numPoints;
    ReshapeSize(3) = [];
    ReshapeSize  = ReshapeSize';
    Spectra = zeros(count,size(XList,2));
    l=1;
    wdf.StartChunkwiseReading();
    while (wdf.AreMoreChunks());
        z = wdf.GetNextDataChunk();
        Size = size(z);
        numberReadIn = Size(1);
        Spectra(l:l+numberReadIn-1,:) = z;%sum(z2);%cumtrapz(z,2);
        l = l + numberReadIn;
    end

a=find(XList>wavenumber-2&XList<wavenumber+2);a=a(1);

ChemicalMap2D=Spectra(:,a);
ChemicalMap3D=reshape(ChemicalMap2D,ReshapeSize(1),ReshapeSize(2),1);

end