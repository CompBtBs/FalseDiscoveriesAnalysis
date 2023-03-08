%
%	sampling_chrr.m  Sampling ENGRO1 and ENGRO2 models by using the CHRR algorithm
%
%	Version 1.0 February 2023
%
%	Authors: 
%		- Bruno G. Galuzzi <bruno.galuzzi@unimib.it> (Department of Biotechnology and Biosciences, University of Milano-Bicocca)
%		- Luca Milazzo <l.milazzo1@campus.unimib.it> (Department of Informatics, Systems, and Communications, University of Milano-Bicocca)
%		- Chiara Damiani <chiara.damiani@unimib.it> (Department of Biotechnology and Biosciences, University of Milano-Bicocca)
%
%	Prerequisites and parameters:
%   	- Install the Cobra Toolbox as described at https://opencobra.github.io/cobratoolbox/stable/installation.html
%		- Fill Line 16 with the path to the installation folder of the Cobra Toolbox (parameter 1)
%       - Set as current directory the folder in which this file is located

addpath('C:\Users\LM856702\Documents\cobratoolbox\') % parameter 1
initCobraToolbox;
modelNames = {'ENGRO 1', 'ENGRO 2'}; 
thinnings = {1, 10, 100};
options.toRound=1;
for modelNamesIndex=1:length(modelNames)
	modelName = modelNames{modelNamesIndex};
    cd('../../samples/')
	mkdir(modelName);
    cd('../models/');
    model = readCbModel(modelName,'fileType','SBML');
    reactions = model.rxns
	model.c=0*model.c;
    cd('../samples/')
    cd(modelName);
	for thinningIndex=1:length(thinnings)
		mkdir(strcat('CHRRThinning',num2str(thinnings{thinningIndex})));
		cd(strcat('CHRRThinning',num2str(thinnings{thinningIndex})))
		for nSample=1000:1000:30000
            rownames = compose('%d', 0:nSample-1);
			options.nStepsPerPoint=thinnings{thinningIndex};
			options.nPointsReturned=nSample;
			for h=0:19
				[modelsampling,samples] = sampleCbModel(model,[], 'CHRR',options); 
				runtime=toc;
				filename=strcat(pwd(),'\',num2str(nSample),'_',num2str(h),'_chrr' ,'.csv');
                table = array2table(samples.', 'VariableNames', reactions, 'RowNames', rownames);
				writetable(table, filename,'WriteRowNames',true);
			end
		end
		cd("..");
    end
    cd("..")
end