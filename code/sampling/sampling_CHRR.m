%
%	sampling_chrr.m  Sampling ENGRO1 and ENGRO2 models by using the CHRR algorithm
%
%   Version 1.0 February 2023
%
%	Authors: 
%		- Bruno Giovanni Galuzzi<> ()
%		- Luca Milazzo <l.milazzo1@campus.unimib.it> (DISCO)
%		- Chiara Damiani <>
%
%	Prerequisites and parameters:
%   	- Install the Cobra Toolbox as described at https://opencobra.github.io/cobratoolbox/stable/installation.html
%		- Fill Line 9 with the path to the installation folder of the Cobra Toolbox (parameter 1)
%		- FIll Line 16 with the path to your sampling folder (parameter 2)
%   

addpath('PATH_TO_COBRA') % parameter 1
initCobraToolbox
modelNames = {'ENGRO1', 'ENGRO2'}; 
thinnings = {1, 10, 100};
options.toRound=1;
for modelNamesIndex=1:length(modelNames)
	modelName = modelNames{modelNamesIndex}
	mydir = strcat('PATH_TO_SAMPLES_FOLDER/' , modelName); % parameter 2
	cd(mydir)
	if(modelName == "ENGRO1"){
		model = readCbModel(strcat('../models/',modelName),'fileType','JSON');
	}else{
		model = readCbModel(strcat('../models/',modelName),'fileType','SMLB');
		listReactions={'EX_O2','EX_Gln','EX_Glc','EX_Arg','EX_THF','EX_Met'};
		model=changeRxnBounds(model,listReactions,[-38,-40,-10,-20,-20,-20],'b');
		model=changeRxnBounds(model,listReactions,[1000],'u');
	}
	model.c=0*model.c;
	for thinningIndex=1:length(thinnings)
		mkdir(strcat("CHRRThinning",num2str(thinnings{thinningIndex})));
		cd(strcat("CHRRThinning",num2str(thinnings{thinningIndex})))
		for nSample=1000:1000:30000
			options.nStepsPerPoint=thinnings{thinningIndex};
			options.nPointsReturned=nSample;
			for h=0:19
				%Generating samples
				[modelsampling,samples] = sampleCbModel(model,[], 'CHRR',options); 
				runtime=toc;
				%Saving the samples
				filename=strcat(pwd(),'\',num2str(nSample),'_',num2str(h),'_chrr' ,'.mat');
				save(filename,'samples');
			end
		end
		cd(mydir)
	end
end