inpath='plmDCAin/';
outpath='plmDCAout/';
interact_path='interact/';
namelist=dir([inpath,'*.txt']);
L=length(namelist);

for i=1:L;
    inpath='plmDCAin/';
    outpath='plmDCAout/';
    interact_path='interact/';
    namelist=dir([inpath,'*.txt']);
    L=length(namelist);
    
    filename=[namelist(i).name];
    name0=filename;
    name1=strsplit(name0,'.');
    
    rna_chain_msa=name1{1};
    splits=strsplit(rna_chain_msa,'_');
    rna=splits{1}
    chain=splits{2};
    rna_chain=strcat(rna,'_',chain);
    
    outname=strcat(outpath, rna_chain_msa,'_evo','.txt');
    
    
    infile=strcat(inpath,filename);
    inname=rna_chain;
    interact_name=strcat(interact_path, rna_chain, '.txt');

    for ii = 10:10:50
        DI_matrix = calculate_evolutionary_constraints_plmDCA_RNAversion(infile, inname, interact_name, ii);
        matrix_name = ['DI_matrix', num2str(ii)];
        DI_Structure.(matrix_name) = sum(DI_matrix, 2);
        DI_column(:,ii/10) = sum(DI_matrix, 2);
    end
    
    writematrix(DI_column, outname,'Delimiter','tab');
    clear
end