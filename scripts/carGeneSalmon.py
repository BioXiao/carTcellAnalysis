import os, sys, re, time, shutil, copy
import pandas as pd
from pprint import pprint

def get_input_list(input_dir):
    input_libs = [ os.path.join(dp, f) \
                   for dp,dn,fn in os.walk(input_dir) \
                   for f in fn if re.search('.*trimmed.fastq$', f) ]
    return input_libs

def get_lib_id(input_lib):
    lib_id = re.search("lib[0-9]+", input_lib).group()
    return lib_id

def get_proj_id(input_lib):
    proj_re = re.compile('P+[0-9]+(-[0-9]+){,1}')
    proj_id = proj_re.search(input_lib).group()

    return proj_id

def get_fc_id(input_lib):
    input_lib = re.sub('EXTERNAL_[A-B]', 'EXTERNAL_', input_lib)
    fc_re = re.compile('((?<=(EXTERNAL_))|(?<=(_[A-B]))).*?XX')
    fc_id = fc_re.search(input_lib)
    if fc_id:
        return fc_id.group()
    else:
        fc_id = re.search('(C|H).*XX', input_lib)
        return fc_id.group()

def run_salmon(input_lib, high_conf=False):
    salmon_exec = "tools/SalmonBeta-0.5.0_OSX-10.10/bin/salmon quant"
    salmon_input = "-r %s" % input_lib

    lib_id = get_lib_id(input_lib)
    tmp_dir = "tmp_%s" % lib_id
    salmon_output = "-o %s" % tmp_dir

    if not high_conf:
        salmon_params = "-i data/indexes/salmon/carTranscript -l U"
    else:
        salmon_params = "-i data/indexes/salmon/GRCh38_CAR -l U"

    salmon_cmd = (' ').join([salmon_exec, salmon_params, salmon_input,
                             salmon_output])
    status = os.system(salmon_cmd)

    if status:
        return "Error!"
    else:
        quant_file = os.path.join(tmp_dir, 'quant.sf')
        quant_df = pd.read_csv(quant_file, sep = '\t', skiprows = 8)
        num_hits = quant_df['NumReads'][0]
        if high_conf:
            tpm = quant_df['TPM'][0]
        else:
            tpm = "NA"
        shutil.rmtree(tmp_dir)
        return (lib_id, num_hits, tpm)

def build_car_gene_dict(input_libs):
    car_gene_dict = {}

    # first pass
    for idx,lib in enumerate(input_libs):
        print "First pass >> lib %s of %s:" % (idx+1, len(input_libs))
        print "Running on %s...\n" % lib
        start_time = time.time()
        lib_id,num_hits,tpm = run_salmon(lib)
        car_gene_dict[lib_id] = {'project': get_proj_id(lib),
                                 'flowcell': get_fc_id(lib),
                                 'num_reads': num_hits, 'tpm': tpm}
        print "%0.2fs elapsed\n" % (time.time() - start_time)

    nonzero_libs = [ l for l in car_gene_dict \
                     if not car_gene_dict[l]['num_reads'] == 0 ]
    input_libs = [ i for i in input_libs \
                   if get_lib_id(i) in nonzero_libs ]

    # second pass
    car_gene_dict_high_conf = copy.deepcopy(car_gene_dict)

    for idx,lib in enumerate(input_libs):
        print "Second pass >> lib %s of %s:" % (idx+1, len(input_libs))
        print "Running on %s...\n" % lib
        start_time = time.time()
        lib_id,num_hits,tpm = run_salmon(lib, high_conf=True)
        car_gene_dict_high_conf[lib_id]['num_reads'] = num_hits
        car_gene_dict_high_conf[lib_id]['tpm'] = tpm
        print "%02.fs elapsed\n" % (time.time() - start_time)

    return (car_gene_dict, car_gene_dict_high_conf)

def write_gene_dict(gene_dict, out_file):
    gene_df = pd.DataFrame(gene_dict).transpose()
    gene_df.to_csv(out_file)

def main(argv):
    input_dir = sys.argv[1]
    output_base = sys.argv[2]

    input_libs = get_input_list(input_dir)
    car_gene_dict,car_gene_dict_high_conf = build_car_gene_dict(input_libs)
    write_gene_dict(car_gene_dict, output_base + '.csv')
    write_gene_dict(car_gene_dict_high_conf, output_base + '_high_conf.csv')


if __name__ == "__main__":
    main(sys.argv[1:])
