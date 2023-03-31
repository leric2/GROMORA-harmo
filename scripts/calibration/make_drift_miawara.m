function [drift, logFile] = make_drift_miawara(rawSpectra,logFile,calibrationTool)

initialIndices_HC={
        find(logFile.isHot & logFile.isMirrorDisp1);      % Hot dis1 
        find(logFile.isHot & logFile.isMirrorDisp2);      % Hot dis2
        find(logFile.isColdSky & logFile.isMirrorDisp1);     % Cold dis1
        find(logFile.isColdSky & logFile.isMirrorDisp2);     % Cold dis2
        };

for i=1:length(initialIndices_HC)
    initialIndices_HC{i}(initialIndices_HC{i}>size(rawSpectra,1))=[];%removal of any indices which dont match raw spectra index
    N_hc(i) = length(initialIndices_HC{i});%number of hot/cold spectra
end

initialIndices_LiRef={
        find(logFile.isLine & logFile.isMirrorDisp1);      % Line        
        find(logFile.isLine & logFile.isMirrorDisp2);      % Line
        find(logFile.isRef & logFile.isMirrorDisp1);     % Reference
        find(logFile.isRef & logFile.isMirrorDisp2);     % Reference
        };

for i=1:length(initialIndices_LiRef)
    initialIndices_LiRef{i}(initialIndices_LiRef{i}>size(rawSpectra,1))=[];%removal of any indices which dont match raw spectra index
    N_lr(i) = length(initialIndices_LiRef{i});%number of sky/ref
end

drift=struct();

for i=1:length(initialIndices_HC);initialIndices_HC{i}=initialIndices_HC{i}(1:min(N_hc)); end; 
for i=1:length(initialIndices_LiRef);initialIndices_LiRef{i}=initialIndices_LiRef{i}(1:min(N_lr)); end; 

drift.dailyMedianHotSpectra1=median(rawSpectra(initialIndices_HC{1},:),1,'omitnan');
drift.dailyStdHotSpectra1=std(rawSpectra(initialIndices_HC{1},:), 'omitnan');

drift.dailyMedianHotSpectra2=median(rawSpectra(initialIndices_HC{2},:),1,'omitnan');
drift.dailyStdHotSpectra2=std(rawSpectra(initialIndices_HC{2},:), 'omitnan');

drift.dailyMedianColdSpectra1=median(rawSpectra(initialIndices_HC{3},:),1,'omitnan');
drift.dailyStdColdSpectra1=std(rawSpectra(initialIndices_HC{3},:),'omitnan');

drift.dailyMedianColdSpectra2=median(rawSpectra(initialIndices_HC{4},:),1,'omitnan');
drift.dailyStdColdSpectra2=std(rawSpectra(initialIndices_HC{4},:),'omitnan');

drift.dailyMedianLineSpectra1=median(rawSpectra(initialIndices_LiRef{1},:),1,'omitnan');
drift.dailyStdLineSpectra1=std(rawSpectra(initialIndices_LiRef{1},:),'omitnan');

drift.dailyMedianLineSpectra2=median(rawSpectra(initialIndices_LiRef{2},:),1,'omitnan');
drift.dailyStdLineSpectra2=std(rawSpectra(initialIndices_LiRef{2},:),'omitnan');

drift.dailyMedianRefSpectra1=median(rawSpectra(initialIndices_LiRef{3},:),1,'omitnan');
drift.dailyStdRefSpectra1=std(rawSpectra(initialIndices_LiRef{3},:),'omitnan');

drift.dailyMedianRefSpectra2=median(rawSpectra(initialIndices_LiRef{4},:),1,'omitnan');
drift.dailyStdRefSpectra2=std(rawSpectra(initialIndices_LiRef{4},:),'omitnan');

ih = [ initialIndices_HC{1}, initialIndices_HC{2}];
ic = [ initialIndices_HC{3}, initialIndices_HC{4}];
il = [ initialIndices_LiRef{1}, initialIndices_LiRef{2}];
ir = [ initialIndices_LiRef{3}, initialIndices_LiRef{4}];

%Outliers
%Check angles are within bounds for hot/cold
hotAngleOutlier=reshape((abs(logFile.Elevation_Angle(ih)-calibrationTool.elevationAngleHot) > calibrationTool.elevationAngleHotTol),[],1);
coldAngleOutlier=reshape((abs(logFile.Elevation_Angle(ic)-calibrationTool.elevationAngleCold) > calibrationTool.elevationAngleColdTol),[],1);

%Check that adc overload does not exceed permitted amount
FFT_adc_overload_hot = reshape((logFile.acqiris_ADC_overflow_counts(ih) > calibrationTool.adcOverloadThresh),[],1);
FFT_adc_overload_cold = reshape((logFile.acqiris_ADC_overflow_counts(ic) > calibrationTool.adcOverloadThresh),[],1);  

%Check for measurements outside of NumberOfStdDev away from mean values
medStdDevThreshHot1=abs((rawSpectra(initialIndices_HC{1},:)-drift.dailyMedianHotSpectra1))>calibrationTool.hotSpectraNumberOfStdDev*drift.dailyStdHotSpectra1;
medStdDevThreshHot2=abs((rawSpectra(initialIndices_HC{2},:)-drift.dailyMedianHotSpectra2))>calibrationTool.hotSpectraNumberOfStdDev*drift.dailyStdHotSpectra2;

medStdDevThreshCold1=abs((rawSpectra(initialIndices_HC{3},:)-drift.dailyMedianColdSpectra1))>calibrationTool.coldSpectraNumberOfStdDev*drift.dailyStdColdSpectra1;
medStdDevThreshCold2=abs((rawSpectra(initialIndices_HC{4},:)-drift.dailyMedianColdSpectra2))>calibrationTool.coldSpectraNumberOfStdDev*drift.dailyStdColdSpectra2;

medStdDevThreshLine1=abs((rawSpectra(initialIndices_LiRef{1},:)-drift.dailyMedianLineSpectra1))>calibrationTool.skySpectraNumberOfStdDev*drift.dailyStdLineSpectra1;
medStdDevThreshLine2=abs((rawSpectra(initialIndices_LiRef{2},:)-drift.dailyMedianLineSpectra2))>calibrationTool.skySpectraNumberOfStdDev*drift.dailyStdLineSpectra2;

medStdDevThreshRef1=abs((rawSpectra(initialIndices_LiRef{3},:)-drift.dailyMedianRefSpectra1))>calibrationTool.refSpectraNumberOfStdDev*drift.dailyStdRefSpectra1;
medStdDevThreshRef2=abs((rawSpectra(initialIndices_LiRef{4},:)-drift.dailyMedianRefSpectra2))>calibrationTool.refSpectraNumberOfStdDev*drift.dailyStdRefSpectra2;

%check how many outliers in each spectrum - dont use if over threshNumRawSpectraHot amount
outlierDetectHot1 = reshape(sum(medStdDevThreshHot1,2)>calibrationTool.threshNumRawSpectraHot,[],1);
outlierDetectHot2 = reshape(sum(medStdDevThreshHot2,2)>calibrationTool.threshNumRawSpectraHot,[],1);

outlierDetectCold1 = reshape(sum(medStdDevThreshCold1,2)>calibrationTool.threshNumRawSpectraCold,[],1);
outlierDetectCold2 = reshape(sum(medStdDevThreshCold2,2)>calibrationTool.threshNumRawSpectraCold,[],1);

outlierDetectLine1 = reshape(sum(medStdDevThreshLine1,2)>calibrationTool.threshNumRawSpectraLine,[],1);
outlierDetectLine2 = reshape(sum(medStdDevThreshLine2,2)>calibrationTool.threshNumRawSpectraLine,[],1);

outlierDetectRef1 = reshape(sum(medStdDevThreshRef1,2)>calibrationTool.threshNumRawSpectraRef,[],1);
outlierDetectRef2 = reshape(sum(medStdDevThreshRef2,2)>calibrationTool.threshNumRawSpectraRef,[],1);

switch calibrationTool.outlierDectectionType
    case 'standard'
      outlierHot = ([outlierDetectHot1 ; outlierDetectHot2] | hotAngleOutlier | FFT_adc_overload_hot)';
      outlierCold = ([outlierDetectCold1 ; outlierDetectCold2] | coldAngleOutlier | FFT_adc_overload_cold )';
      outlierLine = ([outlierDetectLine1 ; outlierDetectLine2])';
      outlierRef =  ([outlierDetectRef1 ; outlierDetectRef2])';
    case 'NoAngleRemoval'
      outlierHot = ([outlierDetectHot1 ; outlierDetectHot2] )';
      outlierCold = ([outlierDetectCold1 ; outlierDetectCold2])';
      outlierLine = ([outlierDetectLine1 ; outlierDetectLine2])';
      outlierRef =  ([outlierDetectRef1 ; outlierDetectRef2])';
end
 
drift.OutHotInd = ih(find(outlierHot));
drift.OutColdInd = ic(find(outlierCold));
drift.OutLineInd = il(find(outlierLine));
drift.OutRefInd = ir(find(outlierRef));

logFile.isOutlier = zeros(size(logFile.isHot));
logFile.isOutlier([drift.OutHotInd,drift.OutColdInd, drift.OutLineInd, drift.OutRefInd]) = 1;

%TODO - Add outliers to calibratedspectra so that they can be added to
%output file
end
