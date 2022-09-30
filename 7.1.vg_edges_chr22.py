import os
import subprocess
import logging
import tempfile
import gzip
import pipes

# For chr 22
start = 16050075
end = 51244237   # last variant position for chr22
num_vertice_linear_bc = 35780755
num_alt_vertices = 71228
Total_vertices = 35893981

end1 = end - start +1
print(end1)
index = start

# get the variant position
with open('variant_positions_snps_indels_chr22.txt', 'r') as f:
    # removing the new line characters
    variant_positions = [line.rstrip() for line in f]

with open('chr22_snps_indel_POS_REF_ALT.txt', 'r') as f:
    # removing the new line characters
    vcf_data_list = [line.rstrip().split() for line in f]

"""
    convert vcf list to vcf dic
"""
# Convert list to dictionary
keys = [vcf_data_list[i][0] for i in range(0,len(vcf_data_list))]
vcf_data_dic ={ k:[] for k in keys }

for i, item in enumerate(vcf_data_list):
    k = keys[i]
    vcf_data_dic[k].append((vcf_data_list[i][1], vcf_data_list[i][2]))

# Add edges accociated with the linear backbone
index = start
with open('linear_bc_chr22_in_variant_range.fa', 'r') as f:
   for line in f:
      lines = f.read()  # skip the first line
      # lines = f.readlines()[1:]  # skip the first line
      
label_list =[]
for label in lines:
    label_list.append(label)

with open('chr22_edges.txt', 'w') as f:
    for i in range(0, end1):
        label = label_list[i]
        edge = (index, index+1, label)
        # print(edge)
        index +=1
        f.write(f"{edge}\n")

 
# Add edges accociated with alternate paths
new_v = num_vertice_linear_bc +1

with open('chr22_edges_2.txt', 'w') as f:
    for pos in variant_positions[0:1700]:
        REF = vcf_data_dic[pos][0][0]
        ALT = vcf_data_dic[pos][0][1]

        ref_len= int(len(str(REF)))
        end_POS_bc= ref_len + int(pos)

        # write outputs to file
        if(int(len(str(ALT))) ==1):
            edge = (pos, end_POS_bc, ALT)
            f.write(f"{edge}\n")

        if(int(len(str(ALT))) !=1):
            for element in ALT:
                if(int(len(str(element))) ==1 and element !=','):
                    edge = (pos, end_POS_bc, element)
                    f.write(f"{edge}\n")
                if(int(len(str(element))) !=1):  
                    v = pos  
                    for i, char in enumerate(element):
                        if(i != len(element)-1):
                            edge = (v, new_v, char)
                            v = new_v
                            new_v +=1
                            f.write(f"{edge}\n")
                        else:
                            edge = (new_v, end_POS_bc, char)
                            f.write(f"{edge}\n")       