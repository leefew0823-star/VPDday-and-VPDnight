clear;clc;
[a,r1]=readgeoraster('F:\result_PKU_NDVI.tif');
[b,r1]=readgeoraster('F:\result_MODIS_NDVI.tif');
[c,r1]=readgeoraster('F:\result_MODIS_EVI.tif');
[d,r1]=readgeoraster('F:\result_PRE.tif');
[e,r1]=readgeoraster('F:\result_Surface.tif');
[f,r1]=readgeoraster('F:\result_ET.tif');
[g,r1]=readgeoraster('F:\result_fix_growseason.tif');
[h,r1]=readgeoraster('F:\result_CSIF.tif');

DCI=nan(360,720);
for i=1:360
    for j=1:720
        aa(1,1)=a(i,j);
        aa(1,2)=b(i,j);
        aa(1,3)=c(i,j);
        aa(1,4)=d(i,j);
        aa(1,5)=e(i,j);
        aa(1,6)=f(i,j);
        aa(1,7)=g(i,j);
        aa(1,8)=h(i,j);
        dci=0;
        for ii=1:8
            if isnan(aa(1,ii))==1
                h=0;
            else
                h=1;
                dci=dci+aa(1,ii)/abs(aa(1,ii));
            end
        end
        DCI(i,j)=dci/8;
    end
    disp(i)
end
