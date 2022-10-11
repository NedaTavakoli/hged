import os
import logging
import tempfile
import gzip
import time
import sys

def create_edges(start, end, num_vertice_linear_bc, variant_pos_file, \
    variant_POS_ALT_REF_file, linear_bc_file, id, linear_edges, alt_edges):

    end1 = end - start +1
    with open(variant_pos_file, 'r') as f:
        variant_positions = [line.rstrip() for line in f]

    variant_positions_index = {k: v+1 for v, k in enumerate(variant_positions)}
        
    with open(variant_POS_ALT_REF_file, 'r') as f:
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
    with open(linear_bc_file, 'r') as f:
        for line in f:
            lines = f.read().rstrip()  # skip the first line
                
        label_list_bc =[]
        for label in lines:
            if label !='\n':
                label_list_bc.append(label)
    index = start
    with open(linear_edges, 'w') as f:
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
    with open(alt_edges, 'w') as f:
        for pos in variant_positions:
            REF = vcf_data_dic[pos][0][0]
            ALT = vcf_data_dic[pos][0][1]
            variant_index = variant_positions_index[str(pos)]
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
if __name__ == "__main__":

    # start = 16050075
    # end = 51244237   # last variant position for chr22
    # num_vertice_linear_bc = 35780755
    # variant_pos_file = variant_positions_snps_indels_chr22.txt
    # variant_POS_ALT_REF_file = chr22_snps_indel_POS_REF_ALT.txt
    # linear_bc_file ='linear_bc_chr22.fa'

    # arguments: start, end, number_of_vertices_linear_backbone, Variant_POS_ALT_REF_file, linear_bc_file
    # id, linear_edges, alt_edges
    start = int(sys.argv[1])
    end = int(sys.argv[2])
    num_vertice_linear_bc = int(sys.argv[3])
    variant_pos_file = sys.argv[4]
    variant_POS_ALT_REF_file= sys.argv[5]
    linear_bc_file= sys.argv[6]
    id= int(sys.argv[7])
    linear_edges=sys.argv[8]
    alt_edges=sys.argv[9]
    print('Creating list of edges...')

    create_edges(start, end, num_vertice_linear_bc, variant_pos_file,\
         variant_POS_ALT_REF_file, linear_bc_file, id, linear_edges, alt_edges)



    