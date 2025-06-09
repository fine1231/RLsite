# RLsite
RLsite: Integrating 3D-CNN and BiLSTM for RNA-Ligand Binding Site Prediction

**Script Requirments:**

* [closeness-degree](https://github.com/KailiWang1/RLBind/tree/main/demo/features)
* [plmDCA](http://zhaoserver.com.cn/ZHmol/ZHmolReSTasite/ZHmolReSTasite.html)
* [SASA](http://yanglab.nankai.edu.cn/RNAsol/)

This is a demo of the input features construction protocol with output format description.

* [SpaceEnergy+LN]

INPUT: RNA structure file
OUTPUT: model input features file

STEP 1: prepare RNA fasta file and RNA monomer structure file
Run ./bin/RLsite 1uud.mol2 -c y
INPUT:  - 1uud.mol2           - RNA structure
OUTPUT: - 1uud#.ism           - RNA space 
        - LN.dat              - RNA Laplace Norm
        - 1uud_sn             - list of local spatial nucleotides
STEP 2: 

The order of steps 2 to 4 can be changed or these steps can be done simultaneously
STEP 2: calculate two network topological properties (closeness and degree) by MATLAB (see refs [1] for more details)
1) calculate RNA static contact network files and nodes closeness centrality using the RNA monomer structure by a default cutoff.
Run closeness/main.m
INPUT:  - closeness/Inputs/1uud.pdb             - RNA monomer structure
OUTPUT: - closeness/Outputs/1uud_mapping.txt      - Output file residues vs PDB file residues
      - closeness/Outputs/1uud_contact.dat      - Static contact network by cutoff
      - closeness/Outputs/1uud_closeness.txt     - Nodes closeness centrality
2) calculate the degree of a graph based on the adjacency matrix.
Run degree/main.m
INPUT:  - degree/Inputs/1uud_contact.dat        - Adjacency matrix file
OUTPUT: - degree/Outputs/1uud_degree.txt        - Nodes degree centrality

STEP 3: calculate predicted accessible surface areas
1) Download the RNAsol software or used the online software (https://yanglab.nankai.edu.cn/RNAsol/) to calculate the solvent accessible surface areas of each nucleotide in RNAs.
INPUT:  - fasta/1uud_B.fasta        - RNA fasta files
OUTPUT: - features/ASA/seq_asa.pre  - RNA accessible surface area files

STEP 4: calculate the evolutional conservation scores for each nucleotide in RNA (http://zhaoserver.com.cn/ZHmol/ZHmolReSTasite/ZHmolReSTasite.html)
1) Obtain the homologous sequences similar to the selected RNA structure of a researched sequence using the BLASTN with the E-value of 0.001 by searching against a non-redundant nucleotide database. You can install BLASTN and MARS database.
INPUT:  - fasta/1uud_B.fasta                  - RNA fasta file
OUTPUT: - MSAoutput/1uud_B.txt    - RNA multiple sequence alignment file
2) Calculate the evolutionary conservation scores using plmdca algorithm.
INPUT:  - plmDCAin/1uud_B_msa.txt                   - RNA monomer structure
OUTPUT: - plmDCAout/1uud_B_msa_evo.txt              - RNA evolutional conservation file

STEP 5: RNA input features integration
Run  ./demo/feature_example/combine.sh

STEP 6: Predicting

The predicting can be executed by:
```
python main.py > ./Output/result.dat
```

The predicted result of RLsite will be output

