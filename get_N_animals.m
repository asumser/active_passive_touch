function [N_animals,AnimalNames]=get_N_animals(vpmIDs,pomIDs,RecDB)
%all, joint, vpm, pom
ANvpm=RecDB{vpmIDs,'AnimalName'};
ANpom=RecDB{pomIDs,'AnimalName'};

AnimalNames{1}=unique(cat(1,ANvpm,ANpom));
AnimalNames{2}=unique(intersect(ANvpm,ANpom));
AnimalNames{3}=unique(ANvpm);
AnimalNames{4}=unique(ANpom);
N_animals=cellfun(@numel,AnimalNames);

