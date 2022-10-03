import os
import subprocess
import logging
import tempfile
import gzip
import pipes
import argparse

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

    # Create variant position dictionary
variant_positions_index = {k: v+1 for v, k in enumerate(variant_positions)}


    
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
        lines = f.read().rstrip()  # skip the first line
            #lines = f.readlines()[1:]  # skip the first line
        
            
    label_list_bc =[]
    for label in lines:
        if label !='\n':
            label_list_bc.append(label)


# Add edges: linear backbone
index = start
with open('chr22_edges.txt', 'w') as f:
    for i in range(0, end1):
        label = label_list_bc[i]
        # get the variant index if it exists
        if  variant_positions_index.__contains__(str(index)) ==1:
            variant_index = variant_positions_index[str(index)]
            edge = (str(index), str(index+1), label, str(variant_index))
        else: # position is not a variant position
            edge = (str(index), str(index+1), label, '-')
        index +=1
        f.write(edge[0]+' '+edge[1]+' '+edge[2]+' '+edge[3]+'\n')


# Add edges accociated with alternate paths
new_v = num_vertice_linear_bc +1

with open('chr22_edges_2.txt', 'w') as f:
    for pos in variant_positions:
        REF = vcf_data_dic[pos][0][0]
        ALT = vcf_data_dic[pos][0][1]
        variant_index = variant_positions_index[str(pos)]
        # variant_index = [i+1 for i, t in enumerate(vcf_data_list) if t[0]== pos]
        ref_len= int(len(str(REF)))
        end_POS_bc= ref_len + int(pos)
        # write outputs to file
        if(int(len(str(ALT))) ==1):
            edge = (str(pos), str(end_POS_bc), ALT, str(variant_index))
            f.write(edge[0]+' '+edge[1]+' '+edge[2]+' '+edge[3]+'\n')
        if(int(len(str(ALT))) !=1):
            ALT_elements = ALT.split(",")
            for element in ALT_elements:
                if(int(len(str(element)))==1):
                    edge = (str(pos), str(end_POS_bc), element,str(variant_index))
                    f.write(edge[0]+' '+edge[1]+' '+edge[2]+' '+edge[3]+'\n')
                else:
                    v = pos  
                    for i, char in enumerate(element.rstrip()):
                        if(i != len(element)-1):
                            edge = (str(v), str(new_v), char, str(variant_index))
                            v = new_v
                            new_v +=1
                            f.write(edge[0]+' '+edge[1]+' '+edge[2]+' '+edge[3]+'\n')
                        else:  # the last element
                            edge = (str(new_v), str(end_POS_bc), char, str(variant_index))
                            f.write(edge[0]+' '+edge[1]+' '+edge[2]+' '+edge[3]+'\n')


# if __name__ == "__main__":


#     parser = argparse.ArgumentParser(description='Get the list of edges for the variation graph')
#     parser.add_argument('--chr', type=int, default='22',
#                         help='Type inteher number to show chromosome number, (default: 22)')
#     parser.add_argument('--start', type=int, default='16050075',
#                         help='Type integer number for the first variant positionm in chromosome (default: 16050075)')
#     parser.add_argument('--end', type=int, default='51244237',
#                         help='Type integer number for the last variant positionm in chromosome (default: 16050075)')
#     parser.add_argument('--bc', type=int, default='51244237',
#                         help='Type integer for number of vertices in linear backbone in chromosome (default: 35780755)')
                    
                
#     args, unparsed = parser.parse_known_args()
#     print('err:', args)

#     id=args.chr
#     start=args.start
#     end=args.end
#     num_vertice_linear_bc=args.bc

#     end1 = end - start +1
#     print(end1)
#     index = start
    