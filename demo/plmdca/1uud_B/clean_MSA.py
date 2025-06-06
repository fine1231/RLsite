import os

inpath = 'MSAoutput'
namelist = [ f for f in os.listdir(inpath) if f.endswith('.txt')]
outpath = 'plmDCAin/'

for item in namelist:
    print(item)
    rna = item.split('_')[0]
    chain = item.split('_')[1]
    
    itempath = os.path.join(inpath, item)
    keeps = []
    with open (itempath, 'r') as file:
        content = file.readlines()
        for index, line in enumerate(content):
            if line.startswith('>'):
                nextindex = index + 1
                string = '/1-' + str(len(content[nextindex]))
                if '(' in line:
                    keep = line.split('(')[0]
                elif index == 0:
                    keep = line.split()[0]
                keep1 = str(keep) + string
                keeps.append(keep1)
            else:
                keeps.append('\n'+line)        
         
         
    outkeep = []            
    for index, line in enumerate(keeps):
        if line.startswith('>'):
                nextindex = index + 1
                nameline = line
                contline = keeps[nextindex]
                if 'N' in keeps[nextindex] or 'R' in keeps[nextindex] or 'Y' in keeps[nextindex] or 'S' in keeps[nextindex] or 'W' in keeps[nextindex] or 'K' in content[nextindex] or 'M' in keeps[nextindex]:
                    pass
                else:
                    outkeep.append(nameline)
                    outkeep.append(contline)
    output = os.path.join(outpath, rna + '_' + chain + '_msa.txt')                
    with open (output, 'w') as f:
        for line in outkeep:
            f.write(line)    

