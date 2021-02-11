import h5py
from pybedtools import BedTool
import sys
sys.path.insert(0, 'CODE/')
import modeling_utils
from modeling_utils import *


def extract_binding_TFs(h5_filepath, out_filepath):
    
    with h5py.File(h5_filepath, 'r') as f, open(out_filepath, 'w') as out:
        print('keys: ', list(f.keys()))

        print('tf_bind keys: ', list(f['tf_binding'].keys()))

        for tf in list(f['tf_binding'].keys()):
            out.write('{}\n'.format(tf))


    print("wrote TF with binding locations data to {}".format(out_filepath))


def test(h5_filepath):
    
    with h5py.File(h5_filepath, 'r') as f:#, open(tfs_file, 'r') as tfs:
        print('keys: ',  list(f.keys()))

        print('tf_bind keys: ', list(f['tf_binding'].keys()))

        for dset in f['tf_binding'].keys():

            dset_data = f['tf_binding'][dset]
            print(dset_data.shape, dset_data.dtype)

            print(dset[:5])

    bed = BedTool(tf_bedfile)
    print(bed.features())
    print(bed.head(50))

def xx(h5_filepath):
    
    with h5py.File(h5_filepath, 'r') as f:
        genes = f['genes'][:].tolist()
        print(len(genes))


if __name__ == '__main__':
    file1 = sys.argv[1] 
    file2 = sys.argv[2]

    #print("file1: {}\nfile2: {}".format(file1,file2))

    common_genes, h_map, c_map = create_gene_index_map(file1,file2,True,'gene_ensg')
    print(len(common_genes))
    print(len(h_map))
    print(len(c_map))
    #extract_binding_TFs(file1, file2)

